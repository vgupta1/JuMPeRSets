###
# UCS Oracle 
###
# At the moment only supports cutting planes
#using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform

using Gurobi

export UCSOracle
export suppFcnUCSRd

#To Do
# Check bounds and use closed-form support function
# Incorporate side-constraints
# Implement reformulation
# Write tests

type UCSOracle <: AbstractOracle
    cons::Dict{Int, UncConstraint} #[cnst index master problem => unc. cnst]
    setup_done::Bool

    eps_kappa::Float64
    Gamma1::Float64
    Gamma2::Float64

    # Cutting plane algorithm
    muhat::Vector{Float64}
    covbar::Array{Float64, 2}
    C::Array{Float64, 2}  #C'C = covbar
    cut_tol::Float64  ##defaults to 1e-6
    cut_model::Model
    cut_vars::Vector{Variable}
    unbounded_support::Bool   #True if no support bounds on any uncertainties

    # Other options
    debug_printcut::Bool
end

function UCSOracle(muhat, covhat, Gamma1, Gamma2, eps_; cut_tol = 1e-6)
    covbar = covhat + Gamma2 * eye(size(covhat)...)
    C = chol(covbar, :U)   #C'C = covbar
    UCSOracle(  Dict{Int, UncConstraint}(), false, 
                eps_, Gamma1, Gamma2, muhat, covbar, C, 
                cut_tol, Model(solver=GurobiSolver(OutputFlag=false)), Variable[], false, false)  
end

#preferred interface
function UCSOracle(data, eps_, delta1, delta2; numBoots=10000, cut_tol=1e-6)
    muhat  = vec(mean(data, 1))
    covhat = cov(data)
    Gamma1 = boot_mu(data, delta1, numBoots)
    Gamma2 = boot_sigma(data, delta2, numBoots)
    UCSOracle(muhat, covhat, Gamma1, Gamma2, eps_, cut_tol=cut_tol)
end


kappa(eps_) = sqrt(1./eps_ - 1.)

#the supp fcn when support = Rd
function suppFcnUCSRd(xs, muhat, Gamma1, Gamma2, covbar, eps_k, cut_sense=:Max)
    toggle = 1.
    if cut_sense == :Min
        xs = copy(-xs)
        toggle = -1.
    end
    norm_x = norm(xs)
    sig_term = sqrt(xs' * covbar * xs)[1]
    ustar =  muhat + Gamma1/norm_x * xs 
    ustar += kappa(eps_k) * sig_term / norm_x / norm_x * xs
    dot(ustar, xs)*toggle, ustar
end

# JuMPeR alerting us that we're handling this contraint
function registerConstraint(w::UCSOracle, con, ind::Int, prefs)
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")
    w.cons[ind] = con

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    haskey(prefs, :cut_tol) && ( w.cut_tol = prefs[:cut_tol] )
    return [:Cut    =>  true]
end

function setup(w::UCSOracle, rm::Model)
    w.setup_done && return
    rd = JuMPeR.getRobust(rm)
    w.cut_model.solver   = rd.cutsolver == nothing ? rm.solver : rd.cutsolver
    d = size(w.covbar, 1)
    @assert (rd.numUncs == d) "Num Uncertainties doesn't match columns in data"

    #w.cut_model.colCat   = rd.uncCat  #only supports continuous variables for now
    @assert (length(rd.uncertaintyset) == 0) #does not support additional cnsts on unctertainties for now

    #if there are no bounds, can use the simpler support function
    ##VG To add
    w.unbounded_support = true
    for i = 1:d
        if (rd.uncLower[i] > -Inf) || (rd.uncUpper[i] < Inf )
            w.unbounded_support=false
            break
        end
    end

    #else build an SOCP for separation
    @defVar(w.cut_model, rd.uncLower[i] <= us[i=1:d] <= rd.uncUpper[i])
    @defVar(w.cut_model, z1[1:d])
    @defVar(w.cut_model, z2[1:d])
    for i = 1:d
        setName(us[i], rd.uncNames[i])
        addConstraint(w.cut_model, us[i] == w.muhat[i] + z1[i] + dot(w.C[:, i], z2))
    end

    addConstraint(w.cut_model, dot(z1, z1) <= w.Gamma1^2 )
    addConstraint(w.cut_model, dot(z2, z2) <= kappa(w.eps_kappa)^2)
    w.cut_vars = us[:]
    w.setup_done = true
end

function generateCut(w::UCSOracle, rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # Update the cutting plane problem's objective
    con = w.cons[ind]
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, m.colVal) #m.colVal= master_sol
    
    #Unscramble
    indx, coeffs = zip(unc_obj_coeffs...)
    d = size(w.covbar, 1)
    xs = zeros(length(coeffs))
    xs[[indx...]] = [coeffs...]

    if w.unbounded_support
        zstar, ustar = suppFcnUCSRd(xs, w.muhat, w.Gamma1, w.Gamma2, w.covbar, w.eps_kappa, cut_sense)
    else
        setObjective(w.cut_model, cut_sense, sum([xs[ix] * w.cut_vars[ix] for ix=1:d]))
        cut_solve_status = solve(w.cut_model, suppress_warnings=true)
        cut_solve_status != :Optimal && error("Cutting plane problem infeasible or unbounded!")
        ustar = getValue(w.cut_vars)
        zstar = getObjectiveValue(w.cut_model)
    end
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

#VG shouldn't be called
function generateReform(w::UCSOracle, rm::Model, ind::Int, master::Model)
	return false
end