#### 
# UI Oracle tests
####
include("../ddusets.jl")  #VG should be dropped
using Base.Test, DDUSets
using JuMPeR

include("test_helpers.jl")

function suppFcnTest()
	srand(8675309)
	data = rand(200, 2)
	zstar, ustar = suppFcnUI([1, 1], data, [0, 0], [1, 1], log(1/.1), .001)
	@test_approx_eq_eps(zstar, 1.6436648065224289, 1e-8)
	@test_approx_eq_eps(ustar[1], .8286168721973264, 1e-8)
	@test_approx_eq_eps(ustar[2], .8150479343251025, 1e-8)

	zstar, ustar = suppFcnUI([1, -1], data, [0, 0], [1, 1], log(1/.1), .001)
	@test_approx_eq_eps(zstar, 0.6764338553682907, 1e-8)
	@test_approx_eq_eps(ustar[1], .837635189604830, 1e-8)
	@test_approx_eq_eps(ustar[2], .16120133423653965, 1e-8)

	zstar, ustar = suppFcnUI([-1, 1], data, [0, 0], [1, 1], log(1/.1), .001)
	@test_approx_eq_eps(zstar, 0.6456428391904294, 1e-8)
	@test_approx_eq_eps(ustar[1], .1708990840068779, 1e-8)
	@test_approx_eq_eps(ustar[2], .8165419231973073, 1e-8)
end

function portTest()
	srand(8675309)
	data = rand(500, 2)
	w = UIOracle(data, [0., 0.], [1., 1.], .1, .2) 
	portTest(w, 0.12335160142638069, [0.5639829612402062, 0.4360170387597938])
end

Test.with_handler(Test.default_handler) do
	suppFcnTest()
	portTest()
end
