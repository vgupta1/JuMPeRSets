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
	@test_approx_eq(zstar, -4.724022297688992)
	@test_approx_eq(ustar[1], -4.213785691942581)
	@test_approx_eq(ustar[2], -0.5102366057464114)

	zstar, ustar = suppFcnFB([1, 1], mfs, mbs, sigfs, sigbs, log_eps, :Max)
	@test_approx_eq(zstar, 5.098525912188081)
	@test_approx_eq(ustar[1], 1.0597051824376162)
	@test_approx_eq(ustar[2], 4.038820729750465)

	zstar, ustar = suppFcnFB([-1,  1], mfs, mbs, sigfs, sigbs, log_eps, :Min)
	@test_approx_eq(zstar, -2.7492629560940403)
	@test_approx_eq(ustar[1], 2.0194103648752324)
	@test_approx_eq(ustar[2], -0.7298525912188081)

	zstar, ustar = suppFcnFB([-1,  1], mfs, mbs, sigfs, sigbs, log_eps, :Max)
	@test_approx_eq(zstar, 6.319708517540586)
	@test_approx_eq(ustar[1], -3.0848542587702927)
	@test_approx_eq(ustar[2], 3.234854258770293)
end

function portTest()
	srand(8675309)
	data = randn(500, 2)
	w = FBOracle(data, .1, .2)
	portTest(w, -1.7936707516055366, [0.585907940563643, 0.41409205943635696])
end

Test.with_handler(Test.default_handler) do
	suppFcnTest()
	portTest()
end
