## PortfolioExperiments 

include("../ddusets.jl")  #VG should be dropped

module PE
using DDUSets, JuMPeR, MLBase
import Distributions

#For N = 500, GammaLCX = 0.2775971931202961
function simMktChen(N; seed=8675309)
	seed != nothing && srand(seed)
	Bs = [.5 * (1 + i/11) for i = 1:10] 
	lbs = -sqrt((1-Bs).*Bs)./(1-Bs)
	ubs = sqrt((1-Bs).*Bs)./Bs
	R = zeros(N, 10)

	for i = 1:10
		d = Distributions.Bernoulli(Bs[i])
		R[:, i] = R[:, i] + ((ubs[i] - lbs[i]) .* rand(d, N) + lbs[i])
	end
	R, lbs, ubs
end

function summarizeChenMkt(Bs)
	lbs = -sqrt((1-Bs).*Bs)./(1-Bs)
	ubs = sqrt((1-Bs).*Bs)./Bs
	map(x-> round(x, 2), [1-Bs lbs ubs])
end

function portTest(oracle, lbs, ubs)
	m = RobustModel() #what is the default solver?
	d = length(lbs)	

	@defUnc(m, lbs[i] <= us[i=1:d] <= ubs[i])
	setDefaultOracle!(m, oracle)

	@defVar(m, xs[1:d] >= 0)
	@addConstraint(m, sum(xs) == 1.)
	@defVar(m, t)
	addConstraint(m, sum([us[i] * xs[i] for i =1:d]) >= t)

	@setObjective(m, Max, t)
	status = solveRobust(m, prefer_cuts=true, report=false)
	@assert status == :Optimal "Did not solve"
	getObjectiveValue(m), getValue(xs)
end

# function portTestCC(R)
# 	m = Model()
# 	N, d = size(R)
# 	@defVar(m, xs[1:d] >= 0)
# 	@addConstraint(m, sum{xs[i], i=1:d} == 1)
# 	@defVar(m, t)
# 	for j=1:N
# 		@addConstraint(m, t <= sum{R[j, i] * xs[i], i=1:d})
# 	end

# 	@setObjective(m, Max, t)
# 	status = solve(m)
# 	@assert status == :Optimal "CC did not solve"
# 	getObjectiveValue(m), getValue(xs)
# end

function portTestCC2(R, q)
	m = Model()
	N, d = size(R)
	@defVar(m, xs[1:d] >= 0)
	@addConstraint(m, sum{xs[i], i=1:d} == 1 )
	@defVar(m, t)

	@addConstraint(m, cnstRefs[j=1:N], t-sum{R[j, i]*xs[i], i=1:d}<= 0)
	@setObjective(m, Max, t)
	status = solve(m)
	@assert status== :Optimal "'CC did not solve"

	#iteratively drop the q worst constraints
	numDropped = 0
	while numDropped < q
		duals = map(cref-> getDual(cref), cnstRefs)
		ix = indmax(duals)  #VG Check this
		chgConstrRHS(cnstRefs[ix], 1e6)
		@assert :Optimal == solve(m) "Iteration $numDropped did not solve"
		numDropped = numDropped + 1
	end
	getObjectiveValue(m), getValue(xs)
end

function calcq_loose(N, delta, d, eps_; trace=false)
	f(q) = min(Distributions.cdf(Distributions.Binomial(N, eps_), q+d) * binomial(q+d, d), 1)
	if trace
		for i = 0:iceil(delta*N)
			println("$i \t $(f(i))")
		end
	end

	#find the maximal queue s.t. f(q) <= delta
	#could binary search, but we're lazy and its fast
	q = iceil(delta*N) + 1
	viol =f(q)
	while viol > delta && q > 0
		q = q-1; viol = f(q)
	end
	return q, viol
end

function each_delta(file_path, N, deltas; seed=8675309)
	f = open(file_path, "w")
	R, lbs, ubs = simMktChen(N, seed=seed)
	writecsv(f, ["Run" "Method" "delta" "score"])

	const k = 5  #k-fold cross validation
	const eps_ = .1
	const d = 10

	#used by many of the procedures to assess
	function test(w, inds)
		zstar, xvals = portTest(w, lbs, ubs)
		quantile(R[inds, :]*xvals[:], eps_) 
	end

	function testCC(w, inds_out)
		z, x = portTestCC2(w, q)
		quantile(R[inds_out, :]*xvals[:], eps_)
	end

	function testNoSupp(w, inds)
		zstar, xvals = portTest(w, fill(-1e6, d), fill(1e6, d))
		quantile(R[inds, :]*xvals[:], eps_) 
	end

	for delta in deltas
		###DY
		try
			scores = cross_validate(inds->UDYOracle(R[inds, :], eps_, delta), 
									test, N, Kfold(N, k))
			writecsv(f, [[1:k] fill("DY", k) fill(delta, k) scores])
			flush(f)
		catch e
			println("DY failed with $e")
		end

		# ###CS 
		try
			scores = cross_validate(inds->UCSOracle(R[inds, :], eps_, delta), 
									test, N, Kfold(N, k))
			writecsv(f, [[1:k] fill("CS", k) fill(delta, k) scores])
			flush(f)
		catch e
			println("CS failed with $e")
		end

		### CS without Support
		try
			scores = cross_validate(inds->UCSOracle(R[inds, :], eps_, delta), 
									testNoSupp, N, Kfold(N, k))
			writecsv(f, [[1:k] fill("CSNSP", k) fill(delta, k) scores])
			flush(f)
		catch e
			println("CSNP failed with $e")
		end

		# ###M
		try
			scores = cross_validate(inds->UMOracle(R[inds, :], lbs, ubs, eps_, delta), 
									test, N, Kfold(N, k))
			writecsv(f, [[1:k] fill("M", k) fill(delta, k) scores])
			flush(f)
		catch e
			println("M failed with $e")
		end

		# ###CC
		q, viol = calcq_loose(N, delta, d, eps_)
		function testCC(w, inds_out)
			z, x = portTestCC2(w, q)
			quantile(R[inds_out, :]*x[:], eps_)
		end

		try 
			scores = cross_validate(inds->R[inds, :], #train
									testCC, N, Kfold(N, k))
			writecsv(f, [[1:k] fill("CC", k) fill(viol, k) scores])
			flush(f)
		catch e 
			println("CC failed with $e")
		end

		# ###LCX
		# #A little cludgey to speed-up testing
		wLCX_= LCXOracle(R[1:(k-1)*N/k, :], eps_, delta)
		GammaLCX = wLCX_.Gamma
		println("Computed GammaLCX: \t $GammaLCX")
		try
			scores = cross_validate(inds->LCXOracle(R[inds, :], eps_, delta, 
									Gamma = GammaLCX, max_iter=5000, 
									ab_cut_tol=5e-5), 
									test, N, Kfold(N, k))
			writecsv(f, [[1:k] fill("LCX", k) fill(delta, k) scores])
			flush(f)
		catch e
			println("LCX failed with $e")
		end
	end
end

#VG Remember to add remarks about DY and CS to show it's equivalent
#VG Add a remark for LCX for same thing.
#Update the LCXOracle tests to work correctly.
#MERGE UP!
#Claen up crap before final merge.

## some trials to figure out which market ot use
function trialLCX(mktGen, N; eps_=.1, delta=.1)
	R, lbs, ubs = mktGen(N)
	wUCS = UCSOracle(R, eps_, delta)
	zCS, xCS = portTest(wUCS, lbs, ubs)

	q, viol = calcq_loose(N, delta, 10, eps_)
	zCC, xCC = portTestCC2(R, q)

	t = @time wLCX = LCXOracle(R, .1, .1, numSamples=1000, numBoots=1000, 
								max_iter=5000, ab_cut_tol=5e-5)
	zLCX, xLCX = portTest(wLCX, lbs, ubs)


	println("Q, viol \t $q, \t $viol")
	println("In Samples: \n LCXs \t $zLCX \n CS \t $zCS \n CC \t $zCC")
	println( map(x->round(x, 2), vcat(xLCX[:]', xCS[:]', xCC[:]')))

	#compute out of Sample Stuf too
	Rout, lbs, ubs = mktGen(5000)
	perfsLCX = Rout * xLCX[:]
	perfsCS  = Rout * xCS[:]
	perfsCC  = Rout * xCC[:]

	println("Quantiles: \n $(quantile(perfsLCX, .1)) \t $(quantile(perfsCS, .1)) \t  $(quantile(perfsCC, .1))")
end

function allPorts(N, GammaLCX; seed=8675309)
	### Solve for all the in-sample portfolios and the zINs
	const eps_ = .1
	const d = 10
	const delta = .1
	R, lbs, ubs = simMktChen(N, seed=seed)
	out = hcat(["Method" "zIn"], ["x$i" for i =1:d]')
	# #CS
	# zstar, xvals = portTest(UCSOracle(R, eps_, delta, lbs, ubs)
	# out = [out; "CS" zstar xvals']

	#CSNP
	zstar, xvals = portTest(UCSOracle(R, eps_, delta), fill(-1000, d), fill(1000, d))
	out = [out; "CSNP" zstar xvals[:]']

	#Ms
	zstar, xvals = portTest(UMOracle(R, lbs, ubs, eps_, delta), lbs, ubs)
	out = [out; "M" zstar xvals[:]']

	# #DY
	# zstar, xvals = portTest(UDYOracle(R, eps_, delta), lbs, ubs)
	# out = [out; "DY" zstar xvals[:]']

	#CC
	q, viol = calcq_loose(N, delta, d, eps_)
	zstar, xvals = portTestCC2(R, q)
	out = [out; "CC" zstar xvals[:]']

	#LCX
	zstar, xvals = portTest(LCXOracle(R, eps_, delta, Gamma=GammaLCX, max_iter=5000), 
										lbs, ubs)
	out = [out; "LCX" zstar xvals[:]']
	out
end

function calcPerf(Rout, xvals; delta=.1)
	perfs = Rout * xvals[:]
	[mean(perfs), std(perfs), quantile(perfs, delta) ]
end

function outSamplePorts(file_path, numRuns, N, GammaLCX, seed=8675309)
	srand(seed)
	const eps_ = .1
	const d = 10
	const delta = .1
	f = open(file_path, "w")

	#Gen an outMkt For assessment purposes
	Rout, lbs, ubs = simMktChen(10000, seed=seed)

	writecsv(f, hcat(["iRun" "Method" "zIn"], ["x$i" for i =1:d]', ["Avg" "Std" "VaR"]) )
	for iRun = 1:numRuns
		R, lbs, ubs = simMktChen(500, seed=nothing)

		#CSNP
		zstar, xvals = portTest(UCSOracle(R, eps_, delta), fill(-1000, d), fill(1000, d))
		writecsv(f, [iRun "CSNP" zstar xvals[:]' calcPerf(Rout, xvals)'])

		#Ms
		zstar, xvals = portTest(UMOracle(R, lbs, ubs, eps_, delta), lbs, ubs)
		writecsv(f, [iRun "M" zstar xvals[:]' calcPerf(Rout, xvals)'])

		#CC
		q, viol = calcq_loose(N, delta, d, eps_)
		zstar, xvals = portTestCC2(R, q)
		writecsv(f, [iRun "CC" zstar xvals[:]' calcPerf(Rout, xvals)'])

		#LCX
		zstar, xvals = portTest(LCXOracle(R, eps_, delta, Gamma=GammaLCX, max_iter=5000), 
											lbs, ubs)
		writecsv(f, [iRun "LCX" zstar xvals[:]' calcPerf(Rout, xvals)'])
		flush(f)
	end
end

end #module

#lydia Davis



#1.2 hours to train an LCX, N=120, delta = .2, Gamma = 0.026812739112985012
#111s per solve... 