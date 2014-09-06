####
# Simple Portfolio Allocation test harness
####
include("UCSOracle.jl")
m = RobustModel()

data = randn(100, 3)
@defUnc(m, us[1:3])
w = UCSOracle(data, 3.0, Gamma1 = .05, Gamma2 = .05)
setDefaultOracle!(m, w)

@defVar(m, xs[1:3] >= 0)
@addConstraint(m, sum(xs) == 1.)
@defVar(m, t)
addConstraint(m, sum([us[i] * xs[i] for i =1:3]) >= t)

@setObjective(m, Max, t)
println("Status: ", solveRobust(m, prefer_cuts=true, report=true, active_cuts=true), "\n")

for ix=1:3
	println( getValue(xs[ix]) )
end






