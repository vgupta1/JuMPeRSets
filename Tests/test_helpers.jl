## Test Helpers for Oracles

using JuMPeR
using Gurobi
function portTest(oracle, zstar_val, ustar_vals; unc_lower=nothing, unc_upper=nothing, TOL=1e-7)
	m = RobustModel(solver=GurobiSolver(OutputFlag=0))
	if unc_lower != nothing
		@defUnc(m, unc_lower[i] <= us[i=1:2] <= unc_upper[i])
	else
		@defUnc(m, us[1:2])
	end

	setDefaultOracle!(m, oracle)
	@defVar(m, xs[1:2] >= 0)
	@addConstraint(m, sum(xs) == 1.)
	@defVar(m, t)
	addConstraint(m, sum([us[i] * xs[i] for i =1:2]) >= t)

	@setObjective(m, Max, t)
	@test :Optimal == solveRobust(m, prefer_cuts=true, report=false, active_cuts=true)
	@test_approx_eq_eps(getObjectiveValue(m), zstar_val, TOL)
	@test_approx_eq_eps(getValue(xs[1]), ustar_vals[1], TOL)
	@test_approx_eq_eps(getValue(xs[2]), ustar_vals[2], TOL)
end
