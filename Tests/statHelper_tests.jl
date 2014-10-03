###
# Helper Tests
###
using Base.Test

include("../helpers.jl")

function boot1Test()
	srand(8765309)
	data = randn(1000)
	@test_approx_eq_eps(boot(data, mean, .9, 10000), 0.10809317304982774, 1e-8)
end

function boot2Test()
	srand(8765309)
	data = randn(1000, 3)
	f= x->mean(minimum(x, 1)) # a nonsense function
	@test_approx_eq_eps(boot(data, f, .9, 10000), -2.8777884239388682, 1e-8)
end

function tapproxTest()
	srand(8765309)
	data = randn(1000, 3)
	out = calcMeansT(data, .1)
	@test_approx_eq_eps(out[1], 0.004451891728022329, 1e-8)
	@test_approx_eq_eps(out[2], 0.06592683689569365, 1e-8)

	out = calcMeansT(data, .1, joint=false)
	@test_approx_eq_eps(out[1], 0.011243060792593858, 1e-8)
	@test_approx_eq_eps(out[2], 0.05913566783112214, 1e-8)
end

function logMeanTest()
	srand(8675309); data = rand(10)
	@test_approx_eq_eps(logMeanExp(2., data), 1.1566289640224618, 1e-8)
	@test_approx_eq_eps(logMeanExp(-1., data), -0.4493802878492036, 1e-8)
end

function fwdBackSigsTest()
	srand(8675309); data = randn(100)
	sigfwd, sigback = calcSigsBoot(data, .1, 10000)
	@test_approx_eq_eps(sigfwd, 1.0918195335954186, 1e-10)
	@test_approx_eq_eps(sigback, 1.08588753435207, 1e-10)
end

function KSGammaTest()
	@test_approx_eq(KSGamma(.1, 1000), 0.0385517413380297)
end

function boot_mu_test()
	srand(8675309); data = randn(100)
	mu_up = boot_mu(data, .1, 100)
	@test_approx_eq(mu_up, 0.15990283885940632)
end

function boot_sigma_test()
	srand(8675309); data = randn(100)
	dist = boot_sigma(data, .1, 100)
	@test_approx_eq(dist, 0.196348437082907)
end

function ab_thresh_test()
	srand(8675309); data = randn(100, 3)
	gamma = calc_ab_thresh(data, .2, 100, 100)	
	@test_approx_eq(gamma, 0.1777841615666656)
end


Test.with_handler(Test.default_handler) do
	ab_thresh_test()
	boot1Test()
	boot2Test()
	tapproxTest()
	logMeanTest()
	fwdBackSigsTest()
	KSGammaTest()
	boot_mu_test()
	boot_sigma_test()
end

