###
# Helper Tests
###
using Base.Test

include("../helpers.jl")

function boot1Test()
	srand(8765309)
	data = randn(1000)
	@test_approx_eq(boot(data, mean, .9, 10000), 0.10809317304982774)
end

function boot2Test()
	srand(8765309)
	data = randn(1000, 3)
	f= x->mean(minimum(x, 1)) # a nonsense function
	@test_approx_eq(boot(data, f, .9, 10000), -2.8777884239388682)
end

function tapproxTest()
	srand(8765309)
	data = randn(1000, 3)
	out = calcMeansT(data, .1)
	@test_approx_eq(out[1], 0.004451891728022329)
	@test_approx_eq(out[2], 0.06592683689569365)

	out = calcMeansT(data, .1, joint=false)
	@test_approx_eq(out[1], 0.011243060792593858)
	@test_approx_eq(out[2], 0.05913566783112214)
end

function logMeanTest()
	srand(8675309); data = rand(10)
	@test_approx_eq(logMeanExp(2., data), 1.1566289640224618)
	@test_approx_eq(logMeanExp(-1., data), -0.4493802878492036)
end

function fwdBackSigsTest()
	srand(8675309); data = randn(100)
	sigfwd, sigback = calcSigsBoot(data, .1, 10000)
	@test_approx_eq(sigfwd, 1.0833792098142592)
	@test_approx_eq(sigback, 1.0836954553273133)
end

Test.with_handler(Test.default_handler) do
	boot1Test()
	boot2Test()
	tapproxTest()
	logMeanTest()
	fwdBackSigsTest()
end

