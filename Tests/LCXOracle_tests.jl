#### 
# LCX Oracle tests
####
include("../ddusets.jl")
using Base.Test, DDUSets
using JuMPeR

include("test_helpers.jl")

function suppFcnTest()
	srand(8675309)
	data = randn(100, 2)
	zstar, ustar = suppFcnLCX([1, 1], data, .1, .35, :Min)
	@test_approx_eq(zstar, -3.76360390834329)
	@test_approx_eq(ustar[1], -1.668344114657504)
	@test_approx_eq(ustar[2], -2.0952597936857864)

	zstar, ustar = suppFcnLCX([1, 1], data, .1, .35, :Max)
	@test_approx_eq(zstar, 4.174542588282224)
	@test_approx_eq(ustar[1], 2.034698238036357)
	@test_approx_eq(ustar[2], 2.139844350245866)

	zstar, ustar = suppFcnLCX([-1, 1], data, .1, .35, :Min)
	@test_approx_eq(zstar, -3.3382572598056868)
	@test_approx_eq(ustar[1], 1.5699603104731068)
	@test_approx_eq(ustar[2], -1.76829694933258)

	zstar, ustar = suppFcnLCX([-1, 1], data, .1, .35, :Max)
	@test_approx_eq(zstar, 3.7547438615152613)
	@test_approx_eq(ustar[1], -1.9828608456376906)
	@test_approx_eq(ustar[2], 1.7718830158775705)
end


function portTest()
	srand(8675309)
	data = randn(500, 2)
	w = LCXOracle(data, .1, .2)
	portTest(w, -1.623657518865325, [0.5374070897836245, 0.4625929102163755])
end

Test.with_handler(Test.default_handler) do
	suppFcnTest()
	portTest()
end
