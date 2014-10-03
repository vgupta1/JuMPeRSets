###
# LCX Oracle 
###
# At the moment only supports cutting planes
using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform
import JuMP.UnsetSolver

export LCXOracle
export suppFcnLCX

#VG To Do
# Merge this and the sampling formulation via an appropriate argument to constructor

#returns bool, violated ab_cut
#Solves 4 LPs
#VG take solver from user
function ab_cut(us, vs, z, data, eps_, CASE; solver=GurobiSolver(OutputFlag=0))
    m = Model(solver=solver)
    @defVar(m, a[1:size(data, 2)])
    @defVar(m, b)

    #define the RHS
    N = size(data, 1)
    @defVar(m, t[1:N] >= 0)
    for i = 1:N
        @addConstraint(m, t[i] >= z*(dot(a, vec(data[i, :])) - b)/N)
    end

    #a,b in a box
    @defVar(m, abs_a[1:size(data, 2)]>=0)
    @defVar(m, abs_b >=0)
    @addConstraint(m, abs_b >= b)
    @addConstraint(m, abs_b >= -b)
    for i = 1:2
        @addConstraint(m, abs_a[i] >=  a[i])
        @addConstraint(m, abs_a[i] >= -a[i])
    end
    #Altered from original paper to match the new sampling technique
    @addConstraint(m, sum{abs_a[i], i=1:size(data, 2)} <= 1.)
    
    if CASE == :1  #VG Do you need the 'vec' calls?
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
function suppFcnLCX(xs, data, eps_, Gamma, cut_sense; 
        bound=1e6, TOL=1e-6, MAXITER=100, trace=false, solver=GurobiSolver(OutputFlag=0))
    m = Model(solver=solver)
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
        iter > MAXITER && error("LCX Supp: Max Iterations exceeded, $(obj - Gamma)")
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
suppFcnLCX(xs, w, cut_sense) = suppFcnLCX(xs, w.data, w.eps_, w.Gamma, cut_sense)

########################
type LCXOracle <: AbstractOracle
    data::Array{Float64, 2}
    eps_::Float64
    Gamma::Float64
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

#Preferred Interface
function LCXOracle(data, eps_, delta; cut_tol=1e-6, debug_printcut=false, 
                    numSamples=int(1e4), numBoots=int(1e4)) 
    @assert (0 < eps_ < 1) "Epsilon invalid: $eps_"
    Gamma = calc_ab_thresh(data, delta, numBoots, numSamples)
    LCXOracle(data, eps_, Gamma, cut_tol, debug_printcut)
end

# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::LCXOracle, rm::Model, ind::Int, prefs) = 
    ! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")


function setup(w::LCXOracle, rm::Model, prefs)
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == size(w.data, 2)) "Num Uncertainties $(rd.numUncs) doesn't match columns in data $(size(w.data, 2))"
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"

    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, w.cut_tol)
end

function generateCut(w::LCXOracle, m::Model, rm::Model, inds::Vector{Int}, active=false)
    new_cons = {}
    rd = JuMPeR.getRobust(rm)

    for ind in inds
        con = JuMPeR.get_uncertain_constraint(rm, ind)
        cut_sense, xs, lhs_const = JuMPeR.build_cut_objective(rm, con, m.colVal) #m.colVal= master_sol
        d = size(w.data, 2)
        zstar, ustar = suppFcnLCX(xs, w, cut_sense)
        lhs_of_cut = zstar + lhs_const

        # SUBJECT TO CHANGE: active cut detection
        if active
            push!(rd.activecuts[ind], 
                JuMPeR.cut_to_scen(ustar, 
                    JuMPeR.check_cut_status(con, lhs_of_cut, w.cut_tol) == :Active))
            continue
        end

        # Check violation
        if JuMPeR.check_cut_status(con, lhs_of_cut, w.cut_tol) != :Violate
            w.debug_printcut && JuMPeR.debug_printcut(rm ,m,w,lhs_of_cut,con,nothing)
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(m, con, ustar)
        w.debug_printcut && JuMPeR.debug_printcut(rm, m, w, lhs_of_cut, con, new_con)
        push!(new_cons, new_con)
    end
    return new_cons
end

#Shouldn't be called
generateReform(w::LCXOracle, m::Model, rm::Model, inds::Vector{Int}) = 0
