###
# UI Oracle 
###
# Only supports cutting planes
using JuMPeR, Optim
import JuMPeR: registerConstraint, setup, generateCut, generateReform

# To Do
# Implement the correct threshold computation
# Write tests

#returns bool, and ustar for degen case
#data is assumed to have two extra rows representing bounds
is_degen(d, Gamma, log_eps) = d * log(1/Gamma) <= log_eps
degen_case(xs, lbounds, ubounds) = [xs[i] >= 0 ? ubounds[i] : lbounds[i]]::Vector{Float64}


#returns zstar, ustar
function suppFcnUI(xs, data_sort, lbounds, ubounds, log_eps, qL, qR; lam_min=1e-8, lam_max = 1e2, xtol=1e-8)
    N = size(data_sort, 1)
    d = size(data_sort, 2)
    Gamma = qL[1]

    if is_degen(d, Gamma, log_eps)
        ustar = degen_case(xs, lbounds, ubounds)
        return dot(xs, ustar), ustar
    end

    #extend the data with the bounds
    const data = [ lbounds' ; data_sort ; ubounds' ]

    function f(lam::Float64)
        term = lam * log_eps
        term2 = 0.0
        for i =1:d
            #attempt to correct for potential overlow
            if xs[i] > 0
                t = exp(xs[i] * (data[:, i] - ubounds[i]) / lam)
                t = dot(qR, t)
                term2 += lam * (xs[i] * ubounds[i]/lam + log(t))
            else
                t = exp(xs[i] * (data[:, i] - lbounds[i]) / lam)
                t = dot(qL, t)
                term2 += lam * (xs[i] * lbounds[i]/lam + log(t))
            end
            term += term2
        end
        term
    end
    res = optimize(f, lam_min, lam_max)
    !res.converged && error("Lambda linesearch did not converge")
    lamstar = res.minimum
    obj1 = res.f_minimum
    ustar = zeros(Float64, d)
    for i = 1:d
        if xs[i] >= 0
            qstar = qR .* exp(xs[i]*data[:, i]/lamstar)
        else
            qstar = qL .* exp(xs[i]*data[:, i]/lamstar)
        end
        qstar /= sum(qstar)
        ustar[i] = dot(qstar, data[:, i])
    end
    dot(ustar, xs), ustar
end

type UIOracle <: AbstractOracle
    cons::Dict{Int, UncConstraint} #[cnst index master problem => unc. cnst]
    setup_done::Bool

    lbounds::Vector{Float64}
    ubounds::Vector{Float64}
    data_sort::Matrix{Float64}
    log_eps::Float64

    # Cutting plane algorithm
    qL::Vector{Float64}
    qR::Vector{Float64}
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

DKWApprox(delta, N) = sqrt( log(2/delta)/N )
function UIOracle(data, lbounds, ubounds, eps, delta; cut_tol = 1e-6) 
    data_sort = zeros(Float64, size(data))
    N, d = size(data)
    for i = 1:d
        data_sort[:, i] = sort(data[:, i])
    end
    Gamma = DKWApprox(delta, N)
    qL = zeros(Float64, N+2)
    qL[1] = Gamma
    qL[2:ifloor(N * (1-Gamma)) + 1] = 1/N
    qL[ifloor(N * (1-Gamma)) + 2] = 1-sum(qL)
    @assert (abs(sum(qL)-1) <= 1e-10) "QL not normalized $(sum(qL))"
    UIOracle(Dict{Int, UncConstraint}(), false, vec(lbounds), vec(ubounds), data_sort, log(1/eps), qL, qL[N+2:-1:1], 1e-6, false)
end

# JuMPeR alerting us that we're handling this contraint
function registerConstraint(w::UIOracle, con, ind::Int, prefs)
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")
    w.cons[ind] = con

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    haskey(prefs, :cut_tol) && ( w.cut_tol = prefs[:cut_tol] )
    return [:Cut    =>  true]
end

function setup(w::UIOracle , rm::Model)
    w.setup_done && return
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == size(w.data_sort, 2)) "Num Uncertainties doesn't match columns in data"

    #ignore any additional constraints on uncertainties for now
    #w.cut_model.colCat   = rd.uncCat
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"
    w.setup_done = true
end

function generateCut(w::UIOracle   , rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # Update the cutting plane problem's objective
    con = w.cons[ind]
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, m.colVal) #m.colVal= master_sol
    
    #Unscramble
    indx, coeffs = zip(unc_obj_coeffs...)
    d = size(w.data_sort, 2)
    xs = zeros(length(coeffs))
    xs[[indx...]] = [coeffs...]

    zstar, ustar = suppFcnUI(xs, w.data_sort, w.lbounds, w.ubounds, w.log_eps, w.qL, w.qR)
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
function generateReform(w::UIOracle, rm::Model, ind::Int, master::Model)
	return false
end