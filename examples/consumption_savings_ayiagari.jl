
using NoLib

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
    # r = 1.025
    # y = 1.0 # available income
    
    K = 20.0
    α = 0.36
    A = 1
    δ = 0.025
    r = 1+α*(1/K)^(1-α) - δ
    w = (1-α)*K^α

    y = w

    c = 0.9*y


    e = 0
    cbar = c


    m = SLVector(;w,r,e)
    s = SLVector(;y)
    x = SLVector(;c)
    y = SLVector(;K)
    z = SLVector(;z=0.0)

    p = SLVector(;β, γ, σ, ρ, cbar, α, δ)

    mc = rouwenhorst(3,ρ,σ)
    
    P = mc.p
    ## decide whether this should be matrix or smatrix
    Q = [SVector(w,r,e) for e in mc.state_values] 

    N = 50

    grid = SGrid(Q) × CGrid(((0.01,4.0,N),))
    
    name = Val(:ayiagari)

    (;name, m, s, x, y, z, p, P, grid)


end


function transition(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, p)
    y = exp(M.e)*M.w + (s.y-x.c)*M.r
    return SLVector( (;y) )
end


function arbitrage(mod::typeof(model), m::SLArray, s::SLArray, x::SLArray, M::SLArray, S::SLArray, X::SLArray, p)
    lhs = 1 - p.β*( X.c/x.c )^(-p.γ)*M.r # >= 0
    rhs = s.y - x.c # >= 0
    res = lhs ⟂ rhs
    return SLVector( (;res) )
end



(;m,s,x,p) = model
M = m
S = s
X = x

transition(model, m, s, x, M, p)
arbitrage(model, m, s, x, M, S, X, p)


sol = NoLib.time_iteration_3(model; verbose=false, improve=false)

using Plots


x0 = sol.solution

P = NoLib.transition_matrix(model, sol.solution)

μ = NoLib.ergodic_distribution(model, sol.solution)



using ForwardDiff
using FiniteDiff



# using Plots
# xvec = [e[1] for e in model.grid[1,:]]
# plot(xvec,μ[1,:])
# plot!(xvec,μ[2,:])
# plot!(xvec,μ[3,:])
# # how to plot distribution ?

# ###

# model.p.r = 1.025

# sol_2 = NoLib.time_iteration_3(model; verbose=false, improve=false)
# μ_2 = NoLib.ergodic_distribution(model, sol_2.solution)
# scatter(xvec,μ_2[1,:])
# scatter!(xvec,μ_2[2,:])
# scatter!(xvec,μ_2[3,:])
# # how to plot distribution ?

# plot([e[2] for e in model.grid[:]],μ[:])


# # sol = NoLib.time_iteration_3(model; verbose=false, improve=true)

# ###

# # Now compute the aggregate equilibrium condition


# s0_l = [NoLib.LVectorLike(merge(model.m, model.s),e)  for e in model.grid[:]]
# x0_l = NoLib.label_GArray(model.x, sol.solution)


# y = LVector(K=40)

# equilibrium(model, sol.solution, μ, y)
# projection(model, y)

# using NoLib: cover




# residual(model, model.y)

# resid(u::Float64) = residual(model, SLVector(K=u))

# Kvec = range(15, 20;length=20)
# Rvec = [resid(e) for e in Kvec]

# using Plots
# plot(Kvec, Rvec)
