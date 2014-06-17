###
# FB Oracle 
###
# At the moment only supports cutting planes
using JuMPeR, Distributions, Optim
import JuMPeR: registerConstraint, setup, generateCut, generateReform
include("bootstrap.jl")

#To Do
# Implement the bootstrapping computation
# Write tests
# Implement reformulation as a SOCP

#Safer versions of log-sum-exp
#VG Update the UI calc to use these
function logSumExp(x, data, probs)
    b = x > 0 ? maximum(data) : minimum(data)
    x*b + log(dot(probs, exp(x * (data-b))))
end

function logSumExp(x, data)
    b = x > 0 ? maximum(data) : minimum(data)
    x*b + log(mean(exp(x * (data-b))))
end

#calculates the means via t approx
# if joint, bounds hold (jointly) simultaneously at level 1-delta_
# o.w. bounds hold individually at level 1-delta_
function calcMeansT(data, delta_; joint=true)
    const N   = length(data)
    const sig = std(data)
    dist      = TDist(N-1)
    delta = joint ? delta_/2 : delta_
    mean(data) + quantile(dist, delta) * sig / sqrt(N), mean(data) + quantile(dist, 1-delta)*sig/sqrt(N)
end

#calculates eq. 21 for a single bootstrap sample
#VG tune the inner optimization for performance
function calcSigSample(boot_sample, CASE)
    mu = mean(boot_sample)
    sig = std(boot_sample)
    b = CASE == :Fwd ? maximum(boot_sample) : minimum(boot_sample)
    f(x) = 2mu/x - 2/x^2 * logSumExp(x, boot_sample)  #include a negative bc we minimize

    if CASE == :Fwd
        res = optimize(f, sig, 10sig)
    elseif CASE == :Back
        res = optimize(f, -10sig, -sig)
    end
    !res.converged && error("Bootstrapping Opt did not converge")

    @assert res.f_minimum < 0
    return sqrt(-res.f_minimum), res.minimum
end

#uses newton method with a hint.
#If we fail, calls other method with warning
function calcSigSample(boot_sample, CASE::Symbol, hint::Float64)
    mu = mean(boot_sample)
    b = CASE == :Fwd ? maximum(boot_sample) : minimum(boot_sample)
    #array version
    f(x) = 2mu/x[1] - 2/x[1]^2 * logSumExp(x[1], boot_sample)  #include a negative bc we minimize

    if CASE == :Fwd
        res = optimize(f, [hint], method=:cg, autodiff=true)
    elseif CASE == :Back
        res = optimize(f, [hint], method=:cg, autodiff=true)
    end
    if !res.f_converged
        warn("Sigsample did not converge with hint.  Running bisection.")
        return calcSigSample(boot_sample, CASE)
    end 

    @assert res.f_minimum < 0
    return sqrt(-res.f_minimum), res.minimum
end

type FwdBackRes
    sigfwd::Float64
    sigback::Float64
    hintfwd::Union(Float64, Nothing)
    hintback::Union(Float64, Nothing)
end

# if joint, bounds hold (jointly) simultaneously at level 1-delta_
# o.w. bounds hold individually at level 1-delta_
# function calcSigsBoot(data::Vector{Float64}, delta_, numBoots; 
#                       CASE=:Both, joint=CASE==:Both, hintfwd = nothing, hintback=nothing)
#     delta = joint ? delta_/2 : delta_
#     sigfwd = 0.; sigback = 0.;
#     if CASE == :Fwd || CASE == :Both
#         if hintfwd == nothing
#             #solve it once on data to get a reasonable starting point
#             dummy, hintfwd = calcSigSample(data, :Fwd)
#         end
#         fwd(boot_sample) = calcSigSample(boot_sample, :Fwd)[1]
#         # fwd(boot_sample) = calcSigSample(boot_sample, :Fwd, hintfwd)[1]
#         sigfwd = boot(data, fwd, 1-delta, numBoots)
#     end
#     if CASE == :Back || CASE == :Both
#         if hintback == nothing
#             dummy, hintback = calcSigSample(data, :Back)
#         end
#         back(boot_sample) = calcSigSample(boot_sample, :Back, hintback)[1]
#         sigback = boot(data, back, 1-delta, numBoots)
#     end
#     FwdBackRes(sigfwd, sigback, hintfwd, hintback)
# end

function calcSigsBoot(data::Vector{Float64}, delta_, numBoots; 
                      CASE=:Both, joint=CASE==:Both)
    delta = joint ? delta_/2 : delta_
    sigfwd = 0.; sigback = 0.;
    if CASE == :Fwd || CASE == :Both
        fwd(boot_sample) = calcSigSample(boot_sample, :Fwd)[1]
        sigfwd = boot(data, fwd, 1-delta, numBoots)
    end
    if CASE == :Back || CASE == :Both
        back(boot_sample) = calcSigSample(boot_sample, :Back)[1]
        sigback = boot(data, back, 1-delta, numBoots)
    end
    FwdBackRes(sigfwd, sigback, nothing, nothing)
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

    return FwdBackRes(sigf, sigb, hintfwd, hintback)
end

#log_eps = log(1/eps_)
#returns zstar, ustar
function suppFcnFB(xs, mfs, mbs, sigfs, sigbs, log_eps, cut_sense)
    toggle = 1
    if cut_sense == "Min"
        xs = copy(-xs)
        toggle = -1
    end
    y1 = zeros(Float64, length(xs))
    y2 = zeros(Float64, length(xs))
    y3 = zeros(Float64, length(xs))
    lam = 0
    for i = 1:length(xs)
        if xs[i] >= 0
            y1[i] = xs[i] * mfs[i]
            y2[i] = xs[i] * sigfs[i]^2
            y3[i] = 0.
            lam += sigfs[i]^2 * xs[i]^2
        else
            y1[i] = xs[i] * mbs[i]
            y2[i] = 0.
            y3[i] = -xs[i] * sigbs[i]^2
            lam += sigbs[i]^2 * xs[i]^2
        end
    end
    lam /= (2 * log_eps)
    lam = sqrt(lam)
    y2 /= lam
    y3 /= lam
    ustar = y1 + y2 - y3
    zstar = dot(xs, ustar) * toggle
    return zstar, ustar
end

type FBOracle <: AbstractOracle
    cons::Dict{Int, UncConstraint} #[cnst index master problem => unc. cnst]
    setup_done::Bool

    mfs::Vector{Float64}
    mbs::Vector{Float64}
    sigfs::Vector{Float64}
    sigbs::Vector{Float64}
    log_eps::Float64

    # Cutting plane algorithm
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

FBOracle(mfs, mbs, sigfs, sigbs, eps; TOL=1e-6) = 
    FBOracle(Dict{Int, UncConstraint}(), true, mfs, mbs, sigfs, sigbs, log(1/eps), TOL, false)

#Preferred constructors
function FBOracle(data, eps, delta1, delta2; CUT_TOL=1e-6, numBoots=1e4)
    N, d  = size(data)
    mfs   = zeros(Float64, d)
    mbs   = zeros(Float64, d)
    sigfs = zeros(Float64, d)
    sigbs = zeros(Float64, d)

    for i = 1:d
        mfs[i], mbs[i] = calcMeansT(data[:, i], delta1/d)
        res = calcSigsBoot(data[:, i], delta2/d, numBoots)
        sigfs[i] = res.sigfwd
        sigbs[i] = res.sigback
    end
    FBOracle(mfs, mbs, sigfs, sigbs, eps, TOL=CUT_TOL)
end
FBOracle(data, eps, delta; CUT_TOL=1e-6, numBoots=1e4) = 
    FBOracle(data, eps, delta/2, delta/2, CUT_TOL=CUT_TOL, numBoots=numBoots)


# JuMPeR alerting us that we're handling this contraint
function registerConstraint(w::FBOracle, con, ind::Int, prefs)
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")
    w.cons[ind] = con

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    haskey(prefs, :cut_tol) && ( w.cut_tol = prefs[:cut_tol] )
    return [:Cut    =>  true]
end

function setup(w::FBOracle , rm::Model)
    w.setup_done && return
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == size(w.mfs, 2)) "Num Uncertainties doesn't match columns in data"
    @assert (length(mfs) == length(mbs) == length(sigfs) == length(sigbs)) "Lengths of means and fb devs dont match uncertainties"

    #ignore any additional constraints on uncertainties for now
    #w.cut_model.colCat   = rd.uncCat
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"
    w.setup_done = true
end

function generateCut(w::FBOracle   , rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # Update the cutting plane problem's objective
    con = w.cons[ind]
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, m.colVal) #m.colVal= master_sol
    
    #Unscramble
    indx, coeffs = zip(unc_obj_coeffs...)
    d = size(w.mfs, 2)
    xs = zeros(length(coeffs))
    xs[[indx...]] = [coeffs...]

    zstar, ustar = suppFcnFB(xs, w.mfs, w.mbs, w.sigfs, w.sigbs, w.log_eps, cut_sense)
    lhs_of_cut = zstar + lhs_const

    #VG Correct this after Iain changes base JumPer
    # TEMPORARY: active cut detection 
    rd = JuMPeR.getRobust(rm)
    if active && (
       ((JuMPeR.sense(con) == :<=) && (abs(lhs_of_cut - con.ub) <= w.cut_tol)) ||
       ((JuMPeR.sense(con) == :>=) && (abs(lhs_of_cut - con.lb) <= w.cut_tol)))
        # Yup its active
        push!(rd.activecuts, ustar[:])
    end
    
    # Check violation
    if !JuMPeR.is_constraint_violated(con, lhs_of_cut, w.cut_tol)
        w.debug_printcut && debug_printcut(rm, m, w, lhs_of_cut,con, nothing)
        return 0  # No new cut
    end
    
    # Create and add the new constraint
    new_con = JuMPeR.build_certain_constraint(con, ustar[:])  #VG correct after Iain updates
    w.debug_printcut && debug_printcut(rm, m, w, lhs_of_cut, con,new_con)
    cb == nothing ? addConstraint(m, new_con) :
                    addLazyConstraint(cb, new_con)
    return 1
end

#Shouldn't be called
function generateReform(w::FBOracle, rm::Model, ind::Int, master::Model)
	return false
end


