####
# Helpers
###
# Contains statistical and other routines used by data-driven sets
using Distributions, Optim, Roots


###Bootstrapping code
# ideally this should all be moved to some base level function
function boot(data::Vector, fun::Function, prob::Float64, numBoots::Int, f_args...)
	const N = size(data, 1)
	dist = DiscreteUniform(1, N)
	out = zeros(Float64, numBoots)
	indices = [1:N]::Vector{Int}
	for i = 1:numBoots
		rand!(dist, indices)
		out[i] = fun(data[indices], f_args...)
	end
	quantile(out, prob)
end

function boot(data::Matrix, fun::Function, prob::Float64, numBoots::Int, f_args...)
	const N = size(data, 1)
	dist = DiscreteUniform(1, N)
	out = zeros(Float64, numBoots)
	indices = [1:N]::Vector{Int}
	for i = 1:numBoots
		rand!(dist, indices)
		out[i] = fun(data[indices, :], f_args...)
	end
	quantile(out, prob)
end

#calculates the means via t approx
# if joint, bounds hold (jointly) simultaneously at level 1-delta_
# o.w. bounds hold individually at level 1-delta_
function calcMeansT(data, delta_; joint=true)
    const N   = length(data)
    const sig_rt_N = std(data)/sqrt(N)
    dist      = TDist(N-1)
    delta = joint ? delta_/2 : delta_
    mean(data) + quantile(dist, delta)*sig_rt_N, mean(data) + quantile(dist, 1-delta)*sig_rt_N
end

#Safer versions of log-sum-exp
function logMeanExp(x::Float64, data_shift::Vector, b::Float64)
    x*b + log(mean(exp(x * (data_shift))))
end
logMeanExp(x::Float64, data::Vector) = (const b = x > 0 ? maximum(data) : minimum(data); logMeanExp(x, data-b, b))

#overwrites hint
function calcSigSampleHint!(boot_sample::Vector{Float64}, CASE::Symbol, hint::Float64, 
						min_u::Float64, max_u::Float64; factor = 2.)
    const mu = mean(boot_sample)

    if CASE == :Fwd
        f(x) = 2mu/x - 2/x^2 * logMeanExp(x, boot_sample-max_u, max_u)  #include a negative bc we minimize
        res = Optim.optimize(f, hint/factor, factor*hint)
    elseif CASE == :Back
        f(x) = 2mu/x - 2/x^2 * logMeanExp(x, boot_sample-min_u, min_u)  #include a negative bc we minimize
        res = Optim.optimize(f, factor*hint, hint/factor)
    end
    !res.converged && error("Bootstrapping Opt did not converge")
    @assert res.f_minimum < 0
    hint = res.minimum
    return sqrt(-res.f_minimum)
end

######
###This is the preferred method
function calcSigsBoot(data::Vector{Float64}, delta_::Float64, numBoots::Int; 
                      CASE=:Both, joint=CASE==:Both)
    const delta  = joint ? delta_/2 : delta_
    sigfwd = 0.; sigback = 0.;
    const min_u::Float64 = minimum(data)
    const max_u::Float64 = maximum(data)

    #Determine an appropriate hint by calling with large params first
    #parameter choices equiv to searching [1e-10, 9std(data)]
    if CASE == :Fwd || CASE == :Both
        hint = 3e-5*std(data)
        calcSigSampleHint!(data, :Fwd, hint, min_u, max_u, factor=3e5*std(data))
        sigfwd = boot(data, calcSigSampleHint!, 1-delta, numBoots, :Fwd, hint, min_u, max_u)
    end
    if CASE == :Back || CASE == :Both
        hint = -3e-5*std(data)
        calcSigSampleHint!(data, :Back, hint, min_u, max_u, factor=3e5*std(data))
        sigback = boot(data, calcSigSampleHint!, 1-delta, numBoots, :Back, hint, min_u, max_u)
    end
    sigfwd, sigback 
end

#could be better about handling overflow
function calcSigsExact(mu, mgf, xmin=1e-10, xmax=1e2)
    f(x) = 2mu/x + 2/x^2 * log( mgf(x) )
    res = optimize(f, xmin, xmax)
    sigf = sqrt(-res.f_minimum)
    hintfwd = res.minimum

    res = optimize(f, -xmax, -xmin)
    sigb = sqrt(-res.f_minimum)
    hintback = res.minimum

    return sigf, sibg
end

#Currently computed using Stephens Approximation 
#Journaly of Royal Statistical Society 1970
function KSGamma(delta, N) 
       const sqrt_N = sqrt(N)
       num = sqrt(.5 * log(2/delta))
       denom = sqrt_N + .12 + .11/sqrt_N
       num/denom
end

function boot_mu(data, delta, numBoots)
    const muhat = mean(data, 1)
    myfun(data_b) = norm(mean(data_b, 1) - muhat)
    boot(data, myfun, 1-delta, numBoots)
end

function boot_sigma(data, delta, numBoots)
    const covhat = cov(data)
    myfun(data_b) = vecnorm(cov(data_b) - covhat)
    boot(data, myfun, 1-delta, numBoots)
end

function bootDY_mu(data, delta, numBoots)
    muhat = mean(data, 1)
    myfun(data_b) = (mu = mean(data_b, 1); ((mu-muhat) * inv(cov(data)) * (mu-muhat)')[1])
    boot(data, myfun, 1-delta, numBoots)
end 

function bootDY_sigma(data, delta, numBoots)
    mu0  = mean(data, 1)
    sig0 = cov(data)
   
    function myfun(data_b)
        muhat = mean(data_b, 1)
        sighat = cov(data_b)
        mudiff = (mu0 - muhat)'*(mu0-muhat)
        f(g2) = eigmin(g2*sighat - sig0 -mudiff)

        #solve
        try       
            fzero(f, [.01, 500])
        catch e
            show(e); println()
            println("Using Manual Bracketing")
            #do your own bracketing
            if f(1) < 0
                lb = 1; ub = 2
                while f(ub) < 0
                    lb = ub
                    ub = ub * 2
                end
            else
                lb = .5; ub = 1
                while f(lb) > 0
                    ub = lb
                    lb = lb * .5
                end
            end
            try
                fzero(f, [lb, ub])
            catch e2
                println("Manual bracketing failed")
                for i = 10.0.^[-2:5]
                    println(f(i))
                end
            end
        end

   end
   boot(data, myfun, 1-delta, numBoots)
end

########
# LCX Stuff
####

########################
# Bootstrapping computation

#Returns indx of last instance of val, 0 if not found
#Assumes sorted_list, and sort_list[start] >= val
function findlast_sort(val::Float64, sort_list::Vector; TOL::Float64=1e-10, start::Int64=1)
    i::Int64 = 1
    for i = int(start):length(sort_list)
        if abs(val - sort_list[i]) > TOL
            break
        end
    end
    i-1
end

## Makes a single pass through zetas/zetahats to solve
# 1/N * max_b { sum_i [zeta_i - b] ^+ - sum_i[zetahat_i - b]^+ }
function singlepass!(zetas::Vector, zetahats::Vector)
    sort!(zetas)  
    sort!(zetahats)

    vstar::Float64 = mean(zetas) - zetas[1] 
    vb::Float64    = mean(zetahats) - zetas[1]

    Gamma::Float64 = vstar - vb
    const N::Int64 = length(zetas)
    pbar::Float64    = 1.0 
    hat_indx::Int64  = 1
    hat_indx_::Int64 = 0
    
    for k = 2:length(zetas)
        vstar += (zetas[k-1] - zetas[k]) * (N-k+1)/N
        hat_indx = findlast_sort(zetas[k-1], zetahats, start=hat_indx_ + 1)
        pbar -=  (hat_indx - hat_indx_)/N
        hat_indx_ = hat_indx
        vb  += (zetas[k-1] - zetas[k]) * pbar
        Gamma = max(Gamma, vstar - vb)
    end
    Gamma
end

#a::Vector is used as storage between bootstraps
function f2(boot_sample::Matrix, data::Matrix, numSamples::Int, a::Vector)
    Gamma::Float64 = 0.
    for i = 1:numSamples
        randn!(a)
        a /= norm(a)
        Gamma = max(Gamma, singlepass!(data*a, boot_sample*a))
    end
    Gamma
end

#Approximates the threshold by sampling a bunch of abs for each bootstrap rep
function calc_ab_thresh(data::Matrix, delta::Float64, numBoots::Int, numSamples::Int)
    const N = size(data, 1)
    const d = size(data, 2)
    a::Vector{Float64} = zeros(Float64, d)
    boot(data, f2, 1-delta, numBoots, data, numSamples, a)
end




