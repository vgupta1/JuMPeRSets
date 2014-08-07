### FwdBackwd Benchmark
using Distributions
include("../FBOracle.jl")

x_dist = Gamma(2, 2)
numBoot = int(1e4)
N = int(1e2)
x_data = min(rand(x_dist, N), 5*std(x_dist))
calcSigsBoot(x_data, .05, numBoot, CASE=:Fwd)        #one for the money
tic()
@profile println( calcSigsBoot(x_data, .05, numBoot, CASE=:Fwd) )  #two for the show
println("Elapsed Time:\t", toc())
#Profile.print()

#baseline  21.8 s
#precomputing 1.14 s
#With a better hint .80
