### 
# Queueing experiments
###
module QExp

using Distributions, Roots, Iterators
include("../helpers.jl")


#these are awkwardly stated in terms of arrival and service rates/ std. of arrivals and serivces
king1(lam, mu, sig_a, sig_s) = ( rho = lam/mu; rho/(1-rho) * .5 * 1/mu * (sig_a^2 * lam^2 + sig_s^2 * mu^2) )
king2(lam, mu, sig_a, sig_s) = ( rho = lam/mu; lam/2/(1-rho) * (sig_a^2 + sig_s^2) )

#more naturally in terms of means and stdsd
kingProb(mu_a, mu_s, sig_a, sig_s, eps) = king1(1/mu_a, 1/mu_s, sig_a, sig_s) / eps

########
#Non-optimized choice of epsilons
########
function FBBound(m_ba, m_fs, sig_ba, sig_fs, eps, n) 
	long_term = log(n/eps) * (sig_fs^2 + sig_ba^2) / 2 / (m_ba-m_fs)
	if m_fs > m_ba  || long_term/(m_ba - m_fs) > n
		return (m_fs - m_ba)*n + sqrt( 2*log(n/eps)*(sig_fs^2 + sig_ba^2) * n )
	end
	return long_term
end

function UCSBound(mu_a, mu_s, sig_a, sig_s, eps, n, Gamma1, Gamma2) 
	long_term = (Gamma1 + sqrt( (n/eps-1)* (sig_a^2 + sig_s^2 + 2Gamma2) ) )^2 / 4 / (mu_a-mu_s)
	if mu_s > mu_a || long_term/(mu_a - mu_s) > n
		return (mu_s - mu_a)*n + (Gamma1 + sqrt( (n/eps-1)* (sig_a^2 + sig_s^2 + 2Gamma2)))*sqrt(n) 
	end
	return long_term
end

##########
#optimizes choice of epsilons for the bound
##########
function FBBound2(m_ba, m_fs, sig_ba, sig_fs, epsbar, n)
	const mubar  = m_fs - m_ba
	const sigbar = sqrt(sig_ba^2 + sig_fs^2)
	function f(t)
		out = -epsbar
		for i = 1:n-1
			out += exp(-.5 * ((t - mubar * i)/sigbar/sqrt(i))^2)
		end
		out
	end
	const f_0 = f(0)
	const fub = FBBound(m_ba, m_fs, sig_ba, sig_fs, epsbar, n)

	if f_0 > 0  #the typical case
		try
			return fzero(f, [0, fub])
		catch e
			show(e); println()
			println("F0 \t $f_0 \t")
			println("Fub \t $fub")
		end
	else  #may happen if not obviously stable
		try
			return fzero(f, fub/2)
		catch e
			show(e); println()
			return Inf
		end
	end
end


function UCSBound2(mu_a, mu_s, sig_a, sig_s, epsbar, n, Gamma1, Gamma2)
	mubar  = mu_s - mu_a
	sigbar = sqrt(sig_a^2 + sig_s^2 + 2Gamma2)
	function f(t)
		out = -epsbar
		for i = 1:n-1
			inside = (t - mubar * i)/(sigbar * sqrt(i)) - Gamma1 / sigbar 
			out += 1/(1 + inside^2)
		end
		out
	end	
	const f_0 = f(0)
	const fub = UCSBound(mu_a, mu_s, sig_a, sig_s, epsbar, n, Gamma1, Gamma2)

	if f_0 > 0  #the typical case
		return fzero(f, [0, fub])
	else  #may happen if not obviously stable
		try
			return fzero(f, fub/2)
		catch e
			show(e); println()
			return Inf
		end
	end
end


##VG Double check this at some point
function get_bpds(data_a, data_s)
	bpds = Float64[]
	idx = 1
	while idx < length(data_a)
		bpd = 1
		W = max(0, data_s[idx] - data_a[idx])
	    while W > 1e-10 && idx < length(data_a)
	        W = max(0, W + data_s[idx] - data_a[idx])
	        bpd += 1
	        idx += 1
	    end
	    push!(bpds, bpd)
	    idx += 1
	end    
	bpds
end

#adjusted quantiles with DKW adjustment
# With probability 1-delta, we have 
# P(busy_period <= p) <= prob for all (p, prob) in (pts, probs)
function bpd_quants(bpds, delta)
	const Gamma  = KSGamma(delta, length(bpds))
	const thresh = quantile(bpds, 1-Gamma)
	pts = sort(filter(b->b <= thresh, unique(bpds)))
	probs = [mean(bpds .< p) for p in pts] + Gamma
	zip(pts, probs)
end

###
#Assumes that relevant statitics already estimated
function FBBoundBPD(m_ba, m_fs, sig_ba, sig_fs, pts_probs, epsbar, n; TOL=1e-4)
	bnd::Float64 = FBBound2(m_ba, m_fs, sig_ba, sig_fs, epsbar, n)
	#search exhaustively for the minimum
	for (pt, prob) in pts_probs  #want 1 - prob < eps_bar
		if 1-prob >= epsbar - TOL  #only need probabilities suff small
			continue
		end
		if pt > n  #only consider bpd(eps1) <= n
			break
		end
		const bnd_  = FBBound2(m_ba, m_fs, sig_ba, sig_fs, epsbar - 1 + prob, pt)
		bnd = min(bnd_, bnd)
	end
	bnd
end

function FBBoundBPD(data_a, data_s, delta, epsbar, n; 
					split = .5, TOL=1e-4, numBoots=int(1e4))
	const N = length(data_a); const splitN = int(split*N)
	bpds = get_bpds(data_a[1:splitN], data_s[1:splitN])
	pts_probs = bpd_quants(bpds, .5delta)

	#now compute the stats...
	const mu_a  = calcMeansT(data_a[splitN+1:N], .25*.5delta, joint=false)[1]
	const mu_s  = calcMeansT(data_s[splitN+1:N], .25*.5delta, joint=false)[2]
	const sig_a = calcSigsBoot(data_a[splitN+1:N], .25*.5delta, numBoots, CASE=:Back)
	const sig_s = calcSigsBoot(data_s[splitN+1:N], .25*.5delta, numBoots, CASE=:Fwd)
	FBBoundBPD(mu_a, mu_s, sig_a, sig_s, pts_probs, epsbar, n, TOL=TOL)
end

#Assumes that relevant statitics already esimtaed
function UCSBoundBPD(mu_a, mu_s, sig_a, sig_s, 
					 pts_probs, epsbar, n, Gamma1, Gamma2; TOL=1e-4)
	bnd::Float64 = UCSBound2(mu_a, mu_s, sig_a, sig_s, epsbar, n, Gamma1, Gamma2)
	#search exhaustively for the minimum
	for (pt, prob) in pts_probs  #want 1 - prob < eps_bar
		if 1-prob >= epsbar - TOL  #only need probabilities suff small
			continue
		end
		if pt > n  #only consider bpd(eps1) <= n
			break
		end
		const bnd_  = UCSBound2(mu_a, mu_s, sig_a, sig_s, epsbar - 1 + prob, pt, Gamma1, Gamma2)
		bnd = min(bnd_, bnd)
	end
	bnd
end

function UCSBoundBPD(data_a, data_s, delta, epsbar, n; split = .5, TOL=1e-4, numBoots=int(1e4))
	const N = length(data_a); const splitN = int(split*N)
	bpds = get_bpds(data_a[1:splitN], data_s[1:splitN])
	pts_probs = bpd_quants(bpds, .5delta)

	#now compute the stats...
	const mu_a  = mean(dat_a[splitN+1:N])
	const mu_s  = mean(dat_s[splitN+1:N])
	const sig_a  = mean(dat_a[splitN+1:N])
	const sig_s  = mean(dat_s[splitN+1:N])
	const Gamma1 = boot_mu([dat_s[splitN+1:N] dat_a[splitN+1:N]], .5delta, numBoots)
	const Gamma2 = boot_sigma([dat_s[splitN+1:N] dat_a[splitN+1:N]], .5delta, numBoots)
	UCSBoundBPD(mu_a, mu_s, sig_a, sig_s, pts_probs, epsbar, n, Gamma1, Gamma2, TOL=TOL)
end


#simulates the quantiles of a transient waiting time
function simTransWait(data_a, data_s, n, epsilons)
	#generate the waiting times using Lindley
	waits = zeros(Float64, ifloor(length(data_a) / (n-1)))
	for run = 1:length(waits)
		W = 0
		for k = 1:n-1
			W = max(0, W + data_s[(run-1)*(n-1) + k] - data_a[(run-1)*(n-1) + k])
		end
		waits[run] = W
	end
	quantile(waits, 1-epsilons)
end


#Best Parameters For Experiment
#Pareto (1, 1.1)  Exponential(3.05)
#utilization .8939
# 275 obs to decide stability at 90%
#Coef factor: approx 1.4
# Pareto trunc at 15  Exp trunc at 5std
dist_a = Pareto(1, 1.1)
dist_s = Exponential(3.05)
const ubound_a = 15
const ubound_s = 5 * std(dist_s)
const delta = .2

function simArrivals(N)
	dat = rand(dist_a, int(N))
	dat[ dat .> ubound_a ] = ubound_a
	dat
end

function simService(N)
	dat = rand(dist_s, int(N))
	dat[ dat .> ubound_s ] = ubound_s
	dat
end

function exp1(n_grid, epsbar_grid; numBoots = int(1e4), N=1000)
	# Compares median waiting times aross methods as n -> infinity
	#Simulate data and fit the Sets
	dat_a = simArrivals(N); dat_s = simService(N)
	
	const m_ba   = calcMeansT(dat_a, .25*.5delta, joint=false)[1]
	const m_fs   = calcMeansT(dat_s, .25*.5delta, joint=false)[2]
	const sig_ba = calcSigsBoot(dat_a, .25*.5delta, numBoots, CASE=:Back)[2]
	const sig_fs = calcSigsBoot(dat_s, .25*.5delta, numBoots, CASE=:Fwd)[1]

	const mu_a  = mean(dat_a)
	const mu_s  = mean(dat_s)
	const sig_a  = mean(dat_a)
	const sig_s  = mean(dat_s)
	const Gamma1 = boot_mu([dat_s dat_a], .5delta, numBoots)
	const Gamma2 = boot_sigma([dat_s dat_a], .5delta, numBoots)

	#####
	## Used for the truncated bound for UFB
	const splitN = int(N/2)
	bpds = get_bpds(dat_a[1:splitN], dat_s[1:splitN])
	pts_probs = bpd_quants(bpds, .5delta)

	#now compute the stats...
	const mu2_ba  = calcMeansT(dat_a[splitN+1:N], .25*.5delta, joint=false)[1]
	const mu2_fs  = calcMeansT(dat_s[splitN+1:N], .25*.5delta, joint=false)[2]
	const sig2_ba = calcSigsBoot(dat_a[splitN+1:N], .25*.5delta, numBoots, CASE=:Back)[2]
	const sig2_fs = calcSigsBoot(dat_s[splitN+1:N], .25*.5delta, numBoots, CASE=:Fwd)[1]
	##########

	####
	# Used for the truncated boundf or UCS
	######
	#now compute the stats...
	const mu2_a  = mean(dat_a[splitN+1:N])
	const mu2_s  = mean(dat_s[splitN+1:N])
	const sig2_a  = mean(dat_a[splitN+1:N])
	const sig2_s  = mean(dat_s[splitN+1:N])
	const Gamma21 = boot_mu([dat_s[splitN+1:N] dat_a[splitN+1:N]], .5delta, numBoots)
	const Gamma22 = boot_sigma([dat_s[splitN+1:N] dat_a[splitN+1:N]], .5delta, numBoots)


	out = zeros(length(n_grid) * length(epsbar_grid), 9)

	for (ix, (n, epsbar)) in enumerate(product(n_grid, epsbar_grid))
		out[ix, 1] = n
		out[ix, 2] = FBBound(m_ba, m_fs, sig_ba, sig_fs, epsbar, n) 
		out[ix, 3] = UCSBound(mu_a, mu_s, sig_a, sig_s, epsbar, n, Gamma1, Gamma2) 
		out[ix, 4] = FBBound2(m_ba, m_fs, sig_ba, sig_fs, epsbar, n)
		out[ix, 5] = UCSBound2(mu_a, mu_s, sig_a, sig_s, epsbar, n, Gamma1, Gamma2)
		out[ix, 6] = FBBoundBPD(mu2_ba, mu2_fs, sig2_ba, sig2_fs, pts_probs, epsbar, n)
		out[ix, 7] = UCSBoundBPD(mu2_a, mu2_s, sig2_a, sig2_s, pts_probs, epsbar, n, Gamma21, Gamma22)
		out[ix, 8] = kingProb(mu_a, mu_s, sig_a, sig_s, epsbar)
		out[ix, 9] = epsbar
	end
	out
end

function exp1(n_grid, file_path, numRuns; numBoots=int(1e4))
	f = open(file_path, "w")
	writecsv(f, ["iRun" n "FB" "CS" "FB2" "CS2" "FBBPD" "CSBPD" "Kingman" "epsbar"])
	for iRun = 1:numRuns
		writecsv(f, [fill(iRun, length(n_grid)) exp1(n_grid, [.5], numBoots=numBoots)])
		flush(f)
	end

	close(f)
end

####################
# Fix eps = .5, n = 10
# Look at convergence of median time as N increases
function exp2(N_grid; numBoots=1e4)
	out = zeros(length(N_grid), 9)
	for (ix, N) in enumerate(N_grid)
		out[ix, :] = [N exp1([10], [.5], N=N, numBoots=numBoots)]
	end
	out
end

function exp2(N_grid, file_path; numRuns=100)
	srand(8675309)
	f = open(file_path, "w")
	writecsv(f, ["iRun" "N" "n" "FB" "CS" "FB2" "CS2" "FBBPD" "CSBPD" "Kingman" "epsbar"])
	for iRun = 1:numRuns
		writecsv(f, [fill(iRun, length(N_grid)) exp2(N_grid, numBoots=500)])
		flush(f)
	end
end

######################
# Vary epsilon
# n = 10, N = 1000
exp3(eps_grid; numBoots=1e4) = exp1([10], eps_grid, numBoots=numBoots, N=1000)


end #ends module qExp




