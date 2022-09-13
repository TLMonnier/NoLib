
using NoLib
const NL=NoLib

using StaticArrays
using LabelledArrays
using NoLib: SGrid, CGrid, PGrid, GArray
import NoLib: transition, arbitrage
import NoLib: ×, ⟂

using QuantEcon: rouwenhorst

model = let 

    β = 0.96
    γ = 4.0
    σ = 0.1
    ρ = 0.0
    r = 1.02
    w = 1.0
    c = 0.9*w
    
    y = 0
    cbar = c

    λ = 0.0

    m = SLVector(;y)
    s = SLVector(;w)
    x = SLVector(;c λ)
    p = LVector(;β, γ, σ, ρ, r, cbar)

    mc = rouwenhorst(3,ρ,σ)
    
    P = mc.p
    
    ## decide whether this should be matrix or smatrix
    Q = [SVector(e) for e in mc.state_values] 

    N = 50

    grid = SGrid(Q) × CGrid(((0.01,4.0,N),))
    
    name = Val(:consumption_savings)

    (;name, m, s, x, p, P, Q, grid)


end


function transition(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    w = exp(M.y) + (s.w-x.c)*p.r
    return SLVector( (;w) )
end

# ⟂ᶠ(a,b)

function arbitrage(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    eq = 1 - p.β*( X.c/x.c )^(-p.γ)*p.r
    @warn "The euler equation is satisfied only if c<w. If c=w, it can be strictly positive."
    eq2 = λ
    return SLVector( (;res) )
end


## Fix the issue with the complementarity


## Solve using time iteration
sol = NoLib.time_iteration(model; verbose=false, improve=false)

## Plot the decision rule

## Compute the ergodic distribution

μ = NoLib.ergodic_distribution(model, sol.solution)


## Plot the ergodic distribution

xvec = [e[1] for e in model.grid[1,:]]
plot(xvec,μ[1,:])
plot!(xvec,μ[2,:])
plot!(xvec,μ[3,:])


## How does it work?

### Compare $\tau()$ and $\tau_fit$

s = [NL.iti(model.grid)...][10]  # take 10th element of the grid
x = sol.solution[10]     # take 10th element of optimal solution

NL.τ(model, s, x)
NL.τ_fit(model, s, x)

### Create a transition matrix corresponding to the transitions for a given decision rule x0
N = length(model.grid)
P = zeros(N,N)
