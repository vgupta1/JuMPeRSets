### 
# Queueing experiments
###
module QExp

using Distributions, Roots
include("../FBOracle.jl")
include("../UCSOracle.jl")

#these are awkwardly stated in terms of arrival and service rates/ std. of arrivals and serivces
king1(lam, mu, sigT, sigx) = ( rho = lam/mu; rho/(1-rho) * .5 * 1/mu * (sigT^2 * lam^2 + sigx^2 * mu^2) )
king2(lam, mu, sigT, sigX) = ( rho = lam/mu; lam/2/(1-rho) * (sigT^2 + sigX^2) )

#more naturally in terms of means and stdsd
kingProb(mu_a, mu_s, sig_a, sig_s, eps) = king1(1/mu_a, 1/mu_s, sig_a, sig_s) / eps

#Non-optimized choice of epsilons
#TO BE Deprecate
FBBound(mbt, mfx, sigbt, sigfx, eps) = log(1/eps) * (sigfx^2 + sigbt^2) / 2 / (mbt-mfx)
UCSBound(mbt, mfx, sigt, sigx, eps) = (1/eps -1) * (sigt^2 + sigx^2) / 4 / (mbt - mfx)

#optimizes choice of epsilons for the bound
function UCSBound(mbt, mfx, sigt, sigx, epsbar, n)
	mubar  = mfx - mbt
	sigbar = sqrt(sigt^2 + sigx^2)
	function f(t)
		out = -epsbar
		for i = 1:n-1
			out += 1/(1 + ((t-mubar * i)/(sigbar * sqrt(i)))^2)
		end
		out
	end
	maxt = mubar*(n-1) + sqrt((n-1)/epsbar -1) * sigbar * (n-1)
	fzero(f, [0., maxt])
end

#optimizes the choice of epsilons for the bound
function FBBound(mbt, mfx, sigt, sigx, epsbar, n)
	const mubar  = mfx - mbt
	const sigbar = sqrt(sigt^2 + sigx^2)
	function f(t)
		out = -epsbar
		for i = 1:n-1
			out += exp(-.5 * ((t-mubar * i)/sigbar/sqrt(i))^2)
		end
		out
	end
	maxt = mubar * (n-1) + sqrt(2log(n/epsbar)*(n-1))*sigbar
	fzero(f, [0., maxt])
end

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

#Copied from UIOracle... Should be moved to a central location
#sup | Fhat - F | \leq Gamma
DKWApprox(delta, N) = sqrt( log(2/delta)/N )

#adjusted quantiles with DKW adjustment
function bpd_quants(bpds, delta)
	const Gamma  = DKWApprox(delta, length(bpds))
	const thresh = quantile(bpds, 1-Gamma)
	pts = sort(filter(b->b <= thresh, unique(bpds)))
	probs = [mean(bpds .< p) for p in pts] + Gamma
	zip(pts, probs)
end

#Assumes that relevant statitics already esimtaed
function UCSBoundBPD(mu_a, mu_s, sig_a, sig_s, 
					 pts_probs, epsbar, n; TOL=1e-4)
	bnd::Float64 = UCSBound(mu_a, mu_s, sig_a, sig_s, epsbar, n)
	#search exhaustively for the minimum
	for (pt, prob) in pts_probs  #want 1 - prob < eps_bar
		if 1-prob >= epsbar - TOL  #only need probabilities suff small
			continue
		end
		if pt > n  #only consider bpd(eps1) <= n
			break
		end
		const bnd_  = UCSBound(mu_a, mu_s, sig_a, sig_s, epsbar - 1 + prob, pt)
		bnd = min(bnd_, bnd)
	end
	bnd
end

function UCSBoundBPD(data_a, data_s, delta, epsbar, n; split = .5, TOL=1e-4, numBoots=int(1e4))
	const N = length(data_a)
	const splitN = int(split*N)
	bpds = get_bpds(data_a[1:splitN], data_s[1:splitN])
	pts_probs = bpd_quants(bpds, .5delta)

	#now compute the stats...
	const mu_a  = calcMeansT(data_a[splitN+1:N], .25*.5delta, joint=false)[1]
	const mu_s  = calcMeansT(data_s[splitN+1:N], .25*.5delta, joint=false)[2]
	const sig_a = boot_sigma(data_a[splitN+1:N], .25*.5delta, numBoots)
	const sig_s = boot_sigma(data_s[splitN+1:N], .25*.5delta, numBoots)
	UCSBoundBPD(mu_a, mu_s, sig_a, sig_s, pts_probs, epsbar, n, TOL=TOL)
end

###
#Assumes that relevant statitics already estimated
function FBBoundBPD(mu_a, mu_s, sig_a, sig_s, pts_probs, epsbar, n; TOL=1e-4)
	bnd::Float64 = FBBound(mu_a, mu_s, sig_a, sig_s, epsbar, n)
	#search exhaustively for the minimum
	for (pt, prob) in pts_probs  #want 1 - prob < eps_bar
		if 1-prob >= epsbar - TOL  #only need probabilities suff small
			continue
		end
		if pt > n  #only consider bpd(eps1) <= n
			break
		end
		const bnd_  = FBBound(mu_a, mu_s, sig_a, sig_s, epsbar - 1 + prob, pt)
		bnd = min(bnd_, bnd)
	end
	bnd
end

function FBBoundBPD(data_a, data_s, delta, epsbar, n; split = .5, TOL=1e-4, numBoots=int(1e4))
	const N = length(data_a)
	const splitN = int(split*N)
	bpds = get_bpds(data_a[1:splitN], data_s[1:splitN])
	pts_probs = bpd_quants(bpds, .5delta)

	#now compute the stats...
	const mu_a  = calcMeansT(data_a[splitN+1:N], .25*.5delta, joint=false)[1]
	const mu_s  = calcMeansT(data_s[splitN+1:N], .25*.5delta, joint=false)[2]
	const sig_a = calcSigsBoot(data_a[splitN+1:N], .25*.5delta, numBoots, CASE=:Back)
	const sig_s = calcSigsBoot(data_s[splitN+1:N], .25*.5delta, numBoots, CASE=:Fwd)
	FBBoundBPD(mu_a, mu_s, sig_a, sig_s, pts_probs, epsbar, n, TOL=TOL)
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

end #ends module qExp