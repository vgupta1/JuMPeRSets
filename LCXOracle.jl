###
# LCX Oracle 
###
# At the moment only supports cutting planes
using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform
include("bootstrap.jl")

#To Do
# Write tests
# Merge this and the sampling formulation via an appropriate argument to constructor

#returns bool, violated ab_cut
#Solves 4 LPs
#VG take solver from user
function ab_cut(us, vs, z, data, eps_, CASE)
    m = Model(solver=GurobiSolver(OutputFlag=0))
    @defVar(m, a[1:2])
    @defVar(m, b)

    #define the RHS
    N = size(data, 1)
    @defVar(m, t[1:N] >= 0)
    for i = 1:N
        @addConstraint(m, t[i] >= z*(dot(a, vec(data[i, :])) - b)/N)
    end

    #a,b in a box
    @defVar(m, abs_a[1:2]>=0)
    @defVar(m, abs_b >=0)
    @addConstraint(m, abs_b >= b)
    @addConstraint(m, abs_b >= -b)
    for i = 1:2
        @addConstraint(m, abs_a[i] >=  a[i])
        @addConstraint(m, abs_a[i] >= -a[i])
    end
   @addConstraint(m, abs_b + abs_a[1] + abs_a[2] <= 1.)
    
    if CASE == :1
        @addConstraint(m, dot(a, vec(vs)) - (z-1) * b >= 0)
        @addConstraint(m, dot(a, vec(us)) - b >= 0)
        @setObjective(m, :Max, dot(a, vec(vs)) - (z-1)*b + dot(a, vec(us)) - b - sum{t[i], i=1:N})
    elseif CASE == :2
        @addConstraint(m, dot(a, vec(vs)) - (z-1) * b <= 0)
        @addConstraint(m, dot(a, vec(us)) - b >= 0)
        @setObjective(m, :Max, dot(a, vec(us)) - b - sum{t[i], i=1:N})
    elseif CASE == :3
        @addConstraint(m, dot(a, vec(vs)) - (z-1) * b >= 0)
        @addConstraint(m, dot(a, vec(us)) - b <= 0)
        @setObjective(m, :Max, dot(a, vec(vs)) - (z-1)*b - sum{t[i], i=1:N})
    end
    status = solve(m)
    status != :Optimal && error("Case $CASE optimization failed with status $status")
    getObjectiveValue(m), getValue(a), getValue(b)
end

#returns the first violation
function ab_cut(us, vs, z, data, eps_, Gamma, TOL)
    obj, astar, bstar = ab_cut(us, vs, z, data, eps_, :1)
    if obj > Gamma + TOL
        return true, astar[:], bstar
    end

    obj, astar, bstar = ab_cut(us, vs, z, data, eps_, :2)
    if obj > Gamma + TOL
        return true, astar[:], bstar
    end

    obj, astar, bstar = ab_cut(us, vs, z, data, eps_, :3)
    if obj > Gamma + TOL
        return true, astar[:], bstar
    end
    return false, Float64[], 0.0
end

#returns zstar,ustar
#VG Remove the dependency on GurobiSolver() explicitly and take frm user
function suppFcnLCX(xs, data, eps_, Gamma, cut_sense; bound=1e6, TOL=1e-6, MAXITER=100, trace=false)
    m = Model(solver=GurobiSolver(OutputFlag=0))
    @defVar(m, -bound <= us[1:2] <= bound)
    @defVar(m, vs[1:2])
    @defVar(m, 1 <= z <= 1/eps_)
    @setObjective(m, cut_sense, dot(xs, us))
    status = solve(m)
    status != :Optimal && error("LCX support init solve failed: $status")

    uvals = getValue(us)
    vvals = getValue(vs)
    zval  = getValue(z)
    obj, astar, bstar = ab_cut(uvals[:], vvals[:], zval, data, eps_, Gamma, TOL)

    iter::Int64 = 1
    while obj > Gamma + TOL
        iter > MAXITER && error("LCX Supp: Max Iterations exceeded")
        iter += 1

        @defVar(m, t[1:2] >= 0)
        m.internalModelLoaded = false  #VG Eliminate this when time comes

        @addConstraint(m, t[1] >= dot(astar[:], vs) - bstar * z + bstar)
        @addConstraint(m, t[2] >= dot(astar[:], us) - bstar)
        @addConstraint(m, t[1] + t[2] <= mean(max(data * astar[:] - bstar, 0))*z + Gamma)
        
        status = solve(m)
        status != :Optimal && error("LCX support Inner Solver Failed: $status")
        uvals = getValue(us)
        vvals = getValue(vs)
        zval  = getValue(z)
        obj, astar, bstar = ab_cut(uvals[:], vvals[:], zval, data, eps_, Gamma, TOL)

        trace && println("Iter: $iter \t $(obj-Gamma) \t $(astar[:]) \t $(bstar[:])")
    end
    getObjectiveValue(m), getValue(us[:])    
end

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

########################
type LCXOracle <: AbstractOracle
    cons::Dict{Int, UncConstraint} #[cnst index master problem => unc. cnst]
    setup_done::Bool

    data::Array{Float64, 2}
    eps_::Float64
    Gamma::Float64

    # Cutting plane algorithm
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

function LCXOracle(data, eps_; cut_tol=1e-6, debug_printcut=false, 
                    numSamples=int(1e4), numBoots=int(1e4), delta=.2) 
    @assert (0 < eps_ < 1) "Epsilon invalid: $eps_"
    Gamma = calc_ab_thresh(data, delta, numBoots, numSamples)
    LCXOracle(Dict{Int,UncConstraint}(), false, data, eps_, Gamma, cut_tol, debug_printcut)
end

LCXOracle(data, eps_, Gamma; cut_tol=1e-6, debug_printcut=false) = 
    ( @assert (0 < eps_ < 1) "Epsilon invalid: $eps_";
    LCXOracle(Dict{Int,UncConstraint}(), false, data, eps_, Gamma, cut_tol, debug_printcut) )

# JuMPeR alerting us that we're handling this contraint
function registerConstraint(w::LCXOracle, con, ind::Int, prefs)
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")
    w.cons[ind] = con

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    haskey(prefs, :cut_tol) && ( w.cut_tol = prefs[:cut_tol] )
    return [:Cut    =>  true]
end

#only supports continuous variables.  
function setup(w::LCXOracle, rm::Model)
    w.setup_done && return
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == size(w.data, 2)) "Num Uncertainties doesn't match columns in data"

    #ignore any additional constraints on uncertainties for now
    #w.cut_model.colCat   = rd.uncCat
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"

    w.setup_done = true
end

function generateCut(w::LCXOracle, rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # Update the cutting plane problem's objective
    con = w.cons[ind]
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, m.colVal) #m.colVal= master_sol
    
    #Unscramble
    indx, coeffs = zip(unc_obj_coeffs...)
    d = size(w.data, 2)
    xs = zeros(length(coeffs))
    xs[[indx...]] = [coeffs...]

    zstar, ustar = suppFcnLCX(xs, w.data, w.eps_, w.Gamma, cut_sense)
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
function generateReform(w::LCXOracle, rm::Model, ind::Int, master::Model)
	return false
end