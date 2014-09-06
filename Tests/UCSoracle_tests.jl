###
# UCS Oracle Tests
###

include("../ddusets.jl")
using Base.Test, DDUSets
using JuMPeR

include("test_helpers.jl")

function suppFcnTest()
	muhat  = [0.0, 1.2]
	Gamma1 = .002
	Gamma2 = .05
	covbar = [ 1.  0.2 
	           0.2 1.3]

	zstar, ustar = suppFcnUCSRd([1, 1], muhat, Gamma1, Gamma2, covbar, .1, :Max)
	@test_approx_eq(zstar, 6.132331444671241)
	@test_approx_eq(ustar[1], 2.4661657223356204)
	@test_approx_eq(ustar[2], 3.66616572233562)

	zstar, ustar = suppFcnUCSRd([1, 1], muhat, Gamma1, Gamma2, covbar, .1, :Min)
	@test_approx_eq(zstar, -3.7323314446712406)
	@test_approx_eq(ustar[1], -2.4661657223356204)
	@test_approx_eq(ustar[2], -1.266165722335620)
end

#doesn't use any bounds
function portTest()
	srand(8675309)
	data = randn(500, 2)
	w = UCSOracle(data, .1, .1, .1)
	portTest(w, -2.4630938710200345, [0.5253456488147422, 0.4746543511852578])
end


function portTest2()
	srand(8675309)
	data = randn(500, 2)
	w = UCSOracle(data, .1, .1, .1)
	portTest(w, -2.4604223389522866, [0.530985977658841, 0.46901402234115896], TOL=1e-6,
			unc_lower=[-1e6, -1e6], unc_upper=[1e6, 1e6])
end


Test.with_handler(Test.default_handler) do
	suppFcnTest()
	portTest()
	portTest2()
end

