#### 
# FB Oracle tests
####
include("../ddusets.jl")
using Base.Test, DDUSets
using JuMPeR

include("test_helpers.jl")

function suppFcnTest()
	mfs = [.1, .2]
	mbs = [-.05, -.25]
	sigfs = [1., 2.]
	sigbs = [2., .5]
	log_eps = log(1/.1)

	zstar, ustar = suppFcnFB([1, 1], mfs, mbs, sigfs, sigbs, log_eps, :Min)
	@test_approx_eq(zstar, -4.12402229768899)
	@test_approx_eq(ustar[1], -4.1137856919425815)
	@test_approx_eq(ustar[2], -0.010236605746411331)

	zstar, ustar = suppFcnFB([1, 1], mfs, mbs, sigfs, sigbs, log_eps, :Max)
	@test_approx_eq(zstar, 5.098525912188081)
	@test_approx_eq(ustar[1], 1.0597051824376162)
	@test_approx_eq(ustar[2], 4.038820729750465)

	zstar, ustar = suppFcnFB([-1,  1], mfs, mbs, sigfs, sigbs, log_eps, :Min)
	@test_approx_eq(zstar, -2.2492629560940407)
	@test_approx_eq(ustar[1], 2.0194103648752324)
	@test_approx_eq(ustar[2], -.22985259121880813)

	zstar, ustar = suppFcnFB([-1,  1], mfs, mbs, sigfs, sigbs, log_eps, :Max)
	@test_approx_eq(zstar, 6.219708517540586)
	@test_approx_eq(ustar[1], -2.984854258770293)
	@test_approx_eq(ustar[2], 3.234854258770293)
end

function portTest()
	srand(8675309)
	data = randn(500, 2)
	w = FBOracle(data, .1, .2)
	portTest(w, -1.6629346889970285, [0.5271213801091257, 0.4728786198908743])
end

Test.with_handler(Test.default_handler) do
	suppFcnTest()
	portTest()
end
