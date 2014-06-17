###
# Sampled LCX oracle
###
# Builds the LCX Oracle by sampling and then passing off to the general oracle for real work

using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform
include("bootstrap.jl")

#To Do
#Change this to be a wrapper around the general oracle and return the oracle, not u

#builds a sampled version of LCX set using the general oracle
function LCXSet(rm, data, eps_, Gamma, numSamples; seed=nothing)
    seed != nothing && srand(seed)
    d = size(data, 2)
    @defUnc(rm, u[1:d])
    @defUnc(rm, v[1:d])
    @defUnc(rm, 1 <= z <= 1/eps_)
    
    #sample a's uniformly from sphere 
    #use the data to generate b's
    for ix = 1:numSamples
        a = randn(d)
        a /= norm(a)
        bs = data * a
        for b in bs
            @defUnc(rm, t1 >= 0)
            @defUnc(rm, t2 >= 0)
            addConstraint(rm, t1 >= dot(a, v) - b * z + b)
            addConstraint(rm, t2 >= dot(a, u) - b )
            rhs = Gamma + mean(max(data * a - b, 0)) * z
            addConstraint(rm, t1 + t2 <= rhs)
        end
    end
    return u
end
