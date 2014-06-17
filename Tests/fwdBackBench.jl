### FwdBackwd Benchmark
using Distributions
include("../FBOracle.jl")

x_dist = Gamma(2, 2)
numBoot = int(1e4)
N = int(1e4)
x_data = min(rand(x_dist, N), 5*std(x_dist))
calcSigsBoot(x_data, .05, numBoot, CASE=:Fwd)        #one for the money
@time(calcSigsBoot(x_data, .05, numBoot, CASE=:Fwd))  #two for the show
