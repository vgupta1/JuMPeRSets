###
# UCS Oracle 
###
# At the moment only supports cutting planes
using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform
include("bootstrap.jl")

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

    # Other options
    debug_printcut::Bool
end

function boot_mu(data, delta, numBoots)
    const muhat = mean(data, 1)
    myfun(data_b) = norm(mean(data_b, 1) - muhat)
    boot(data, myfun, 1-delta, numBoots)
end

function boot_sigma(data, delta, numBoots)
    const covhat = cov(data)
    myfun(data_b) = normfro(cov(data_b) - covhat)
    boot(data, myfun, 1-delta, numBoots)
end

# data stored in rows
function UCSOracle(data, eps_kappa; Gamma1 = nothing, Gamma2 = nothing, 
                                    delta1 = nothing, delta2 = nothing, 
                                    numBoots= int(1e5), 
                                    cut_tol=1e-6) 
    Gamma1 == nothing && delta1 == nothing && error("Must supply either Gamma1 or delta1")
    Gamma2 == nothing && delta2 == nothing && error("Must supply either Gamma2 or delta2")
    Gamma1 == nothing && (Gamma1 = boot_mu(data, delta1, numBoots))
	Gamma2 == nothing && (Gamma2 = boot_sigma(data, delta2, numBoots))

    println("Gamma1 \t $Gamma1")
    println("Gamma2 \t $Gamma2")

	covbar = cov(data) + Gamma2 * eye(size(data, 2))
    C = chol(covbar, :U)  #C'C = covbar
    UCSOracle(  Dict{Int, UncConstraint}(), false, 
    			eps_kappa, Gamma1, Gamma2, vec(mean(data, 1)), covbar, C, 
                cut_tol, Model(), Variable[], 
                false)  
end

kappa(eps) = sqrt(1/eps - 1)

#the supp fcn when support = Rd
function suppFcnRd(xs, muhat, Gamma1, Gamma2, covbar, eps_k)
    norm_x = norm(xs)
    sig_term = sqrt(xs' * w.covbar * xs)[1]
    ustar =  muhat + Gamma1/norm_x * xs 
    ustar += kappa(eps_k) * sig_term / norm_x / norm_x * xs
    ustar, dot(ustar, xs)
end
suppFcnRd(xs, w::UCSOracle) = suppFcnRd(xs, w.muhat, w.Gamma1, w.Gamma2, w.covbar, w.eps_kappa)

# JuMPeR alerting us that we're handling this contraint
function registerConstraint(w::UCSOracle, con, ind::Int, prefs)
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")
    w.cons[ind] = con

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    haskey(prefs, :cut_tol) && ( w.cut_tol = prefs[:cut_tol] )
    return [:Cut    =>  true]
end

#only supports continuous variables.  
function setup(w::UCSOracle, rm::Model)
    w.setup_done && return
    rd = JuMPeR.getRobust(rm)
    w.cut_model.solver   = rd.cutsolver == nothing ? rm.solver : rd.cutsolver
    d = size(w.covbar, 1)
    @assert (rd.numUncs == d) "Num Uncertainties doesn't match columns in data"


    #ignore any additional constraints on uncertainties for now
    #w.cut_model.colCat   = rd.uncCat
    println( length(rd.uncertaintyset) )
    @assert (length(rd.uncertaintyset) == 0)

    @defVar(w.cut_model, rd.uncLower[i] <= us[i=1:d] <= rd.uncUpper[i])
    @defVar(w.cut_model, z1[1:d])
    @defVar(w.cut_model, z2[1:d])
    for i = 1:d
        setName(us[i], rd.uncNames[i])
        addConstraint(w.cut_model, us[i] == w.muhat[i] + z1[i] + dot(w.C[:, i], z2))
    end

    addConstraint(w.cut_model, dot(z1, z1) <= w.Gamma1 )
    addConstraint(w.cut_model, dot(z2, z2) <= kappa(w.eps_kappa))
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

    setObjective(w.cut_model, cut_sense, sum([xs[ix] * w.cut_vars[ix] for ix=1:d]))
    cut_solve_status = solve(w.cut_model, suppress_warnings=true)
    cut_solve_status != :Optimal && error("Cutting plane problem infeasible or unbounded!")
    ustar = getValue(w.cut_vars)
    lhs_of_cut = getObjectiveValue(w.cut_model) + lhs_const

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