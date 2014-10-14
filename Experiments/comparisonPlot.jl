####
# Supp Fcn Plots
# Just testing for now... to be updated
####
include("../ddusets.jl")  #VG should be dropped
using DDUSets

module SP
using Gurobi, Distributions, JuMPeR

function simData(N, seed = 8675309)
	srand(seed)
	d1 = Beta(4, 4); d2 = Beta(.4, 2)
	[2 * rand(d1, N)-1  2*rand(d2, N)-1]
end

function createSuppFcnPlot(oracle::AbstractOracle; numDirs=200, bounds=nothing)
	out = Array(Float64, numDirs, 5)	
	for ix in 1:numDirs
		m = RobustModel(solver=GurobiSolver(OutputFlag=0), cutsolver=GurobiSolver(OutputFlag=0))
		@defVar(m, t)
		if bounds !=nothing
			lbnds, ubnds = zip(bounds...)
			@defUnc(m, lbnds[i] <= us[i=1:2] <= ubnds[i])
		else
			@defUnc(m, us[i=1:2])
		end
		setDefaultOracle!(m, oracle)

		theta = 2pi / numDirs * ix
		cr = addConstraint(m, cos(theta) * us[1] + sin(theta) * us[2] <= t)
		@setObjective(m, Min, t)
		try
			solveRobust(m, active_cuts=true, prefer_cuts=true, report=false)

			#extract the cut itself
			rd = JuMPeR.getRobust(m)
			ustar = getScenario(cr).data
			out[ix, :] = [cos(theta) sin(theta) getValue(t) ustar']
			println(ix, "\t", out[ix, :])
		catch
			#do nothing
		end
	end
	out
end

function saveExp(cuts, path)
	f = open(path, "w")
	writecsv(f, ["N" "x" "y" "zstar" "u1" "u2"])
	writecsv(f, cuts)
	close(f)
end

#For the LCX plot
# 100 data points, Gamma = .060  (beta densities)
# 1000 data points, GAmma = .017
# 10000 data points, Gamma = 0.005380338502093185... num Boots only 1000

end

dat = SP.simData(10000)
# wLCX = LCXOracle(dat, .1, .2, numBoots=1000)
# wUDY = UDYOracle(dat, .1, .2)
# wUCS = UCSOracle(dat, .1, .2)
