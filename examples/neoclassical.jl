
using NoLib

using StaticArrays
using LabelledArrays
using NoLib: SGrid, CGrid, PGrid, GArray
import NoLib: transition, arbitrage
import NoLib: ×

## Define Model

model = let 

    α = 0.3
    β = 0.96
    γ = 2.0
    δ = 0.1

    k = ((1/β - (1-δ))/α)^(1/(α-1))
    i = δ*k
    z = 0.0

    p = (;α, β, γ, δ)
    
    m = SLVector( (; z))
    s = SLVector( (; k))
    x = SLVector( (; i))

    P = @SMatrix [0.9 0.1; 0.1 0.9]
    # Q = @SMatrix [-0.05; 0.05]
    Q = @SMatrix [-0.1; 0.1]

    exo = SGrid( [Q[i,:] for i=1:size(Q,1)] )
    endo = CGrid( ((0.1, 5.0, 100),) )
    grid = exo × endo


    (;m, s, x, p, P, Q, exo, endo, grid)

end

function transition(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    K = s.k*(1-p.δ) + x.i
    return SLVector( (;K) )
end

function arbitrage(model::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    c = exp(m.z)*s.k^p.α - x.i
    C = exp(M.z)*S.k^p.α - X.i
    r = p.β*(C/c)^(-p.γ)*(1-p.δ + p.α*exp(M.z)*S.k^(p.α-1)) - 1
    return SLVector( (;r) )
end


### Solve the model
res = NoLib.time_iteration(model;improve=false)

### Plot decision rule
φ = res.solution

using Plots
using NoLib: iti

xvec = [e[1] for e in model.grid[2,:]]
yvec = [e[1] for e in φ[2,:]]
plot(xvec, yvec)

# improve the plot....



## Rewriting time_iteration

### Check SGrid, CGrid, PGrid objects

using NoLib: SGrid, CGrid, PGrid, ×

### Check GArray object

### GArray objects represent a vector of points, matching the geometry of the grid
### Create constant initial guess
using Garray
x0 = ...

### Compare On-grid and off-grid indexing (cf interp.jl)
x0[]
x0()

### Check the transition iterator $\tau$

### Write the optimality function `F(model, s, x, x0::GArray)` where s is a grid point, x a controls and  x0 a vector of controls


### Vectorize the optimality function `F(model, x1, x0)`

### Given x1, compute the derivative w.r.t. x0. Solve for the optimal x1 given x0.

### Compute the derivative w.r.t. x1. Make an improvement step.

### Write a time iteration method with an improvement option.


