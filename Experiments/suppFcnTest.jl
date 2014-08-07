####
# Supp Fcn Plots
# Just testing for now... to be updated
####
using Gurobi
include("UCSOracle.jl")
include("LCXSample.jl")
include("LCXOracle.jl")
include("FBOracle.jl")
include("UIOracle.jl")

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
		setDefaultOracle!(m, w)

		theta = 2pi / numDirs * ix
		addConstraint(m, cos(theta) * us[1] + sin(theta) * us[2] <= t)
		@setObjective(m, Min, t)
		solveRobust(m, active_cuts=true, prefer_cuts=true, report=false)

		#extract the cut itself
		rd = JuMPeR.getRobust(m)
		ustar = transpose(rd.activecuts[1][1:2])
		out[ix, :] = [cos(theta) sin(theta) getValue(t) ustar]
	end
	out
end

#unc_factory(rm) = uncs
function createSuppFcnPlot(unc_factory::Function; numDirs=200)
	out = Array(Float64, numDirs, 5)	
	for ix in 1:numDirs
		m = RobustModel(solver=GurobiSolver(OutputFlag=0), cutsolver=GurobiSolver(OutputFlag=0))
		us = unc_factory(m)
		@defVar(m, t)

		theta = 2pi / numDirs * ix
		addConstraint(m, cos(theta) * us[1] + sin(theta) * us[2] <= t)
		@setObjective(m, Min, t)
		solveRobust(m, active_cuts=true, prefer_cuts=true, report=false)

		#extract the cut itself
		rd = JuMPeR.getRobust(m)
		ustar = transpose(rd.activecuts[1][1:2])
		out[ix, :] = [cos(theta) sin(theta) getValue(t) ustar]
	end
	out
end

function testLCX(numSamples, numDirs)
    data = randn(500, 2)
    unc_factoryLCX(model) = LCXSet(model, data, .1, 0.0, numSamples, seed=8675309) 
    out2 = createSuppFcnPlot(unc_factoryLCX, numDirs=numDirs)
end

data = randn(500, 2)
#w = UIOracle(data, [-5., -5.], [5., 5.], .1, .1)
w = FBOracle(data, .1, .2)

out = createSuppFcnPlot(w, numDirs = 20)
println(out)

