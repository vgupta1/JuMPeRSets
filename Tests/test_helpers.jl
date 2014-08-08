## Test Helpers for Oracles

using JuMPeR
function portTest(oracle, zstar_val, ustar_vals)
	m = RobustModel()
	setDefaultOracle!(m, oracle)
	@defUnc(m, us[1:2])

	@defVar(m, xs[1:2] >= 0)
	@addConstraint(m, sum(xs) == 1.)
	@defVar(m, t)
	addConstraint(m, sum([us[i] * xs[i] for i =1:2]) >= t)

	@setObjective(m, Max, t)
	@test :Optimal == solveRobust(m, prefer_cuts=true, report=false, active_cuts=true)
	@test_approx_eq(getObjectiveValue(m), zstar_val)
	@test_approx_eq(getValue(xs[1]), ustar_vals[1])
	@test_approx_eq(getValue(xs[2]), ustar_vals[2])
end
