####
# Supp Fcn Plots
# Just testing for now... to be updated
####
include("../ddusets.jl")  #VG should be dropped
using DDUSets

module SP
using Gurobi, Distributions, JuMPeR

function simData(N)
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
			@defUnc(m, lbbnds[i] <= us[i=1:2] <= ubnds[i])
		else
			@defUnc(m, us[i=1:2])
		end
		setDefaultOracle!(m, oracle)

		theta = 2pi / numDirs * ix
		cr = addConstraint(m, cos(theta) * us[1] + sin(theta) * us[2] <= t)
		@setObjective(m, Min, t)
		solveRobust(m, active_cuts=true, prefer_cuts=true, report=false)

		#extract the cut itself
		rd = JuMPeR.getRobust(m)
		ustar = getScenario(cr).data
		#ustar = transpose(rd.activecuts[1][1:2])
		out[ix, :] = [cos(theta) sin(theta) getValue(t) ustar']
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
# 100 sample points, Gamma = .060  (beta densities)
# 1000 sample points, GAmma = .017


end