###
# FB Oracle 
###
# At the moment only supports cutting planes
import JuMPeR: registerConstraint, setup, generateCut, generateReform

export FBOracle
export suppFcnFB

#log_eps = log(1/eps_)
#returns zstar, ustar
function suppFcnFB(xs, mfs, mbs, sigfs, sigbs, log_eps, cut_sense=:Max)
    sign_flip = 1
    if cut_sense == :Min
        xs = copy(-xs)
        sign_flip = -1
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
    zstar = dot(xs, ustar) * sign_flip
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
function FBOracle(data, eps, delta1, delta2; CUT_TOL=1e-6, numBoots=int(1e4))
    N, d  = size(data)
    mfs   = zeros(Float64, d)
    mbs   = zeros(Float64, d)
    sigfs = zeros(Float64, d)
    sigbs = zeros(Float64, d)

    for i = 1:d
        mfs[i], mbs[i] = calcMeansT(data[:, i], delta1/d)
        sigfs[i], sigbs[i] = calcSigsBoot(data[:, i], delta2/d, numBoots)
    end
    FBOracle(mfs, mbs, sigfs, sigbs, eps, TOL=CUT_TOL)
end
FBOracle(data, eps, delta; CUT_TOL=1e-6, numBoots=int(1e4)) = 
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


