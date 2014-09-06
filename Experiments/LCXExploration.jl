###
#LCXTesting
####
#Some testing/exploration fcns for the LCX Set

module LCXTest  #module space for iterative changes

include("../bootstrap.jl")
include("../LCXOracle.jl")

# #Approximates the threshold by sampling a bunch of abs for each bootstrap rep
# #VG Could be optimized
# function calc_ab_thresh(data::Matrix{Float64}, delta, numBoots, numSamples)
#     const N = size(data, 1)
#     const d = size(data, 2)
#     a = zeros(Float64, d)

#     function boot_func(boot_sample)
#         Gamma::Float64 = 0.0
#         for i = 1:numSamples
#             randn!(a)
#             a /= norm(a)
#             const as = data * a
#             const as_sample = boot_sample * a
#             for b in as
#                 Gamma = max(Gamma, mean(max(as_sample -b, 0)) - mean(max(as-b, 0)))
#             end
#         end
#         Gamma
#     end
#     boot(data, boot_func, 1-delta, numBoots)
# end

# #######################
# ## Some functions for standard gaussian
# function rhs(a, b)
#     norm_a = norm(a)
#     norm_dist = Normal()
#     norm_a * pdf(norm_dist, b/norm_a) - b*cdf(norm_dist, -b/norm_a)
# end

# varisk(sigma, epsilon) = sigma * quantile(Normal(), 1-epsilon)
# cvar(sigma, epsilon) = sigma * pdf(Normal(), quantile(Normal(), 1-epsilon)) / epsilon

# #######################
# function separate(us, vs, z, data, CASE)
#     m = Model(solver=GurobiSolver(OutputFlag=0))
#     @defVar(m, a[1:2])
#     @defVar(m, b)

#     #define the RHS
#     N = size(data, 1)
#     @defVar(m, t[1:N] >= 0)
#     for i = 1:N
#         @addConstraint(m, t[i] >= z*(dot(a, vec(data[i, :])) - b)/N)
#     end

#     #a,b in a box
#     @defVar(m, abs_a[1:2]>=0)
#     @defVar(m, abs_b >=0)
#     @addConstraint(m, abs_b >= b)
#     @addConstraint(m, abs_b >= -b)
#     for i = 1:2
#         @addConstraint(m, abs_a[i] >=  a[i])
#         @addConstraint(m, abs_a[i] >= -a[i])
#     end
#    @addConstraint(m, abs_b + abs_a[1] + abs_a[2] <= 1)
    
#     if CASE == :1
#         @addConstraint(m, dot(a, vec(vs)) - (z-1) * b >= 0)
#         @addConstraint(m, dot(a, vec(us)) - b >= 0)
#         @setObjective(m, :Max, dot(a, vec(vs)) - (z-1)*b + dot(a, vec(us)) - b - sum{t[i], i=1:N})
#     elseif CASE == :2
#         @addConstraint(m, dot(a, vec(vs)) - (z-1) * b <= 0)
#         @addConstraint(m, dot(a, vec(us)) - b >= 0)
#         @setObjective(m, :Max, dot(a, vec(us)) - b - sum{t[i], i=1:N})
#     elseif CASE == :3
#         @addConstraint(m, dot(a, vec(vs)) - (z-1) * b >= 0)
#         @addConstraint(m, dot(a, vec(us)) - b <= 0)
#         @setObjective(m, :Max, dot(a, vec(vs)) - (z-1)*b - sum{t[i], i=1:N})
#     end
#     status = solve(m)
#     status != :Optimal && error("Case $CASE optimization failed with status $status")

#     return (getObjectiveValue(m), getValue(a), getValue(b))
# end

# function separate(us, vs, z, data; TOL=1e-8)
#     #Right now solve all of them an return worst one
#     obj1, astar1, bstar1 = separate(us, vs, z, data, :1)
#     obj2, astar2, bstar2 = separate(us, vs, z, data, :2)
#     obj3, astar3, bstar3 = separate(us, vs, z, data, :3)
    
#     if obj1 >= max(obj1, obj2, obj3) - TOL
#         return obj1, astar1, bstar1
#     elseif obj2 >= max(obj1, obj2, obj3) - TOL
#         return obj2, astar2, bstar2
#     else
#         return obj3, astar3, bstar3
#     end
# end

# #solves via dual simplex
# function suppFcn2(xs, data, Gamma, eps_; bound=1e6, TOL=1e-6, trace=false)
#     m = Model(solver=GurobiSolver(OutputFlag=0))
#     @defVar(m, -bound <= us[1:2] <= bound)
#     @defVar(m, vs[1:2])
#     @defVar(m, 1 <= z <= 1/eps_)

#     @setObjective(m, Max, dot(xs, us))
#     status = solve(m)
#     status != :Optimal && error("Initial Solve Failed")
#     uvals = getValue(us)
#     vvals = getValue(vs)
#     zval = getValue(z)
#     iter = 1
    
#     obj, astar, bstar = separate(uvals[:], vvals[:], zval, data)

#     while obj > Gamma + TOL
#         t1 = max(dot(astar[:], vvals) - bstar * zval + bstar, 0)
#         t2 = max(dot(astar[:], uvals) - bstar, 0)
        
#         @defVar(m, t[1:2] >= 0)
# 		m.internalModelLoaded = false
# 		# MathProgBase.updatemodel!(getInternalModel(m))

#         @addConstraint(m, t[1] >= dot(astar[:], vs) - bstar * z + bstar)
#         @addConstraint(m, t[2] >= dot(astar[:], us) - bstar)
#         @addConstraint(m, t[1] + t[2] <= mean(max(data * astar[:] - bstar, 0))*z + Gamma)

#         status = solve(m)
#         status != :Optimal && error("Inner Solver Failed status $status")
#         uvals = getValue(us)
#         vvals = getValue(vs)
#         zval = getValue(z)
#         obj, astar, bstar = separate(uvals[:], vvals[:], zval, data)

#         trace && println(iter, "  ", obj - Gamma, "  ", astar[:]', "  ", bstar)
#         iter +=1
#     end
#     getObjectiveValue(m), uvals, vvals, zval
# end




##### A sampled version
# function suppFcn(x)
#     m = Model(solver=GurobiSolver())
#     @defVar(m, us[1:2])
#     @defVar(m, vs[1:2])
#     @defVar(m, 1 <= z <= 1/eps_)
#     for ix = 1:numSamples
#         @defVar(m, tu >= 0)
#         @defVar(m, tv >= 0)
#         @addConstraint(m, dot(vec(ab_samples[ix, 1:2]), us) - ab_samples[ix, 3] <= tu)
#         @addConstraint(m, dot(vec(ab_samples[ix, 1:2]), vs) - ab_samples[ix, 3] * (z-1) <= tv)
#         @addConstraint(m, tu + tv <= rhs(ab_samples[ix, 1:2], ab_samples[ix, 3]) * z)
#     end
#    @setObjective(m, Max, dot(x, us))
#     solve(m)

#     #figure out whose tight
#     println("z = $(getValue(z))")
#     uvals = getValue(us)
#     vvals = getValue(vs)
#     zval  = getValue(z)
    
#     tus = max(ab_samples[:, 1:2] * vec(uvals[:]) - ab_samples[:, 3], 0)
#     tvs = max(ab_samples[:, 1:2] * vec(vvals[:]) - ab_samples[:, 3]*(zval -1), 0)
    
#     lhs = [zval * rhs(ab_samples[ix, 1:2], ab_samples[ix, 3]) for ix = 1:numSamples]
#     lhs -= (tus + tvs)

#     for (i, v) in enumerate(lhs)
#         if v < 1e-10
#             println(v, "  ", ab_samples[i, :])
#         end
#     end

#     getObjectiveValue(m), getValue(us)
# end


for logN = 2:5
    #generate some data
    data = randn(int(10^logN), 2)

    for logSam = 2:8
        tic();
        w = LCXOracle(data, .1, numSamples=int(10^logSam), numBoots=int(1e5))
        t = toq();
        println("$logN \t $logSam \t $(w.Gamma) \t $t")
    end
end

end