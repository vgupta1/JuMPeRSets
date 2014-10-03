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
	@test_approx_eq(zstar, -3.1455496412236243)
	@test_approx_eq(ustar[1], -1.313810356077738)
	@test_approx_eq(ustar[2], -1.8317392851458862)

	zstar, ustar = suppFcnLCX([1, 1], data, .1, .35, :Max)
	@test_approx_eq(zstar, 3.447833040263504)
	@test_approx_eq(ustar[1], 1.5387299947906286)
	@test_approx_eq(ustar[2], 1.9091030454728752)

	zstar, ustar = suppFcnLCX([-1, 1], data, .1, .35, :Min)
	@test_approx_eq(zstar, -2.8182398717617616)
	@test_approx_eq(ustar[1], 1.3068007701076567)
	@test_approx_eq(ustar[2], -1.5114391016541049)

	zstar, ustar = suppFcnLCX([-1, 1], data, .1, .35, :Max)
	@test_approx_eq(zstar, 3.157077158917673)
	@test_approx_eq(ustar[1], -1.7212256391067058)
	@test_approx_eq(ustar[2], 1.4358515198109671)
end


function portTest()
	srand(8675309)
	data = randn(500, 2)
	w = LCXOracle(data, .1, .2, 1e-6, false)
	portTest(w, -1.4293015271256384, [0.5424436347803758, 0.4575563652196242])
end

Test.with_handler(Test.default_handler) do
	suppFcnTest()
	portTest()
end
