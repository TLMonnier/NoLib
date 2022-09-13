F(model, s, x::SVector, φ::GArray) = 
    sum(
         w*arbitrage(model,s,x,S,φ(S)) 
         for (w,S) in τ(model, s, x)
    )


F(model, controls::GArray, φ::GArray) =
    GArray(
        model.grid,
        [
            F(model,s,x,φ)
            for (s,x) in zip(iti(model.grid), controls)
        ],
    )

dF_1(model, controls::GArray, φ::GArray) =
    GArray(    # this shouldn't be needed
        model.grid,
        [
            ForwardDiff.jacobian(u->F(model, s, u, φ), x)
            for (s,x) in zip(iti(model.grid), controls) 
        ]
    )

dF_2(model, s, x::SVector, φ::GArray, dφ::GArray) = 
    sum(
            w*ForwardDiff.jacobian(u->arbitrage(model,s,x,S,u), φ(S))* dφ(S)
            for (w, S) in τ(model, s, x)
    )


dF_2(model, controls::GArray, φ::GArray, dφ::GArray) =
    GArray(
        model.grid,
        [(dF_2(model,s,x,φ,dφ)) for (s,x) in zip(iti(model.grid), controls) ],
    )

using LinearMaps

dF_2(model, x::GArray, φ::GArray) = let
    n = length(φ)*length(eltype(φ))
    fun = u->(ravel(dF_2(model, x, φ, unravel(x, u))))
    xxx = ravel(x)
    # return xxx
    e = fun(xxx)
    LinearMap(
        fun,
        n,
        n
    )
end


function time_iteration_1(model;
    T=500,
    K=10,
    tol_ε=1e-8,
    tol_η=1e-6,
    verbose=false
)

    N = length(model.grid)
    x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
    x1 = deepcopy(x0)

    local x0
    local x1
    local t

    for t=1:T

        r0 = F(model, x0, x0)
        ε = norm(r0)

        if ε<tol_ε
            return (;message="Solution found", solution=x0, n_iterations=t)
        end

        if verbose
            println("ϵ=$(ε)")
        end

        x1.data .= x0.data
        
        for k=1:K

            r1 = F(model, x1, x0)
            J = dF_1(model, x1, x0)
            # dF!(J, model, x1, x0)

            dx = GArray(model.grid, J.data .\ r1.data)

            η = norm(dx)
            
            x1.data .-= dx.data

            verbose ? println("Iteration $t: $η") : nothing

            if η<tol_η
                break
            end

        end

        x0 = x1

    end

    return (;solution=x0, message="No Convergence", n_iterations=t)

end

using NLsolve

function time_iteration_2(model;
    T=500,
    K=10,
    tol_ε=1e-8,
    tol_η=1e-6,
    verbose=false
)

    N = length(model.grid)
    x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
    x1 = deepcopy(x0)

    local x0
    local x1
    local t


    for t=1:T

        function fun(u::AbstractVector{Float64})
            x = unravel(x0, u)
            r = F(model, x, x0)
            return ravel(r)
        end

        function dfun(u::AbstractVector{Float64})
            x = unravel(x0, u)
            dr = dF_1(model, x, x0)
            J = convert(Matrix,dr)
            return J
        end

        u0 = ravel(x0)
        sol = nlsolve(fun, dfun, u0)
        u1 = sol.zero

        η = maximum(u->abs(u[1]-u[2]), zip(u0,u1))
        
        verbose ? println("Iteration $t: $η") : nothing
        
        x0 = unravel(x0, u1)

        if η<tol_η
            return (;solution=x0, message="Convergence", n_iterations=t)
        end

    end

    return (;message="No Convergence", n_iterations=T)

end

using LinearAlgebra

function time_iteration(model;
    T=500,
    K=10,
    tol_ε=1e-8,
    tol_η=1e-6,
    verbose=false,
    improve=true
)

    version_check()

    N = length(model.grid)
    x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
    x1 = deepcopy(x0)

    local x0
    local x1
    local t


    for t=1:T

        function fun(u::AbstractVector{Float64})
            x = unravel(x0, u)
            r = F(model, x, x0)
            return ravel(r)
        end

        function dfun(u::AbstractVector{Float64})
            x = unravel(x0, u)
            dr = dF_1(model, x, x0)
            J = convert(Matrix,dr)
            return J
        end

        u0 = ravel(x0)
        sol = nlsolve(fun, dfun, u0)
        u1 = sol.zero

        η = maximum(u->abs(u[1]-u[2]), zip(u0,u1))
        
        verbose ? println("Iteration $t: $η") : nothing
        
        x1 = unravel(x0, u1)

        if !(improve)
            x0 = x1
        else
            # x = T(x)
            # xnn = T(xn)
            # x - xnn = -T'(x) (x - xn)
            # x = xnn - T' (x - xn)
            # x = (I-T')\(xnn - T' xn)

            # TODO: accelerate this part
            A = dF_1(model, x1, x0)
            MA = convert(Matrix, A)
            B = dF_2(model, x1, x0)
            MB = convert(Matrix, B)
            
            Tp = -MA\MB

            u11 = (I-Tp)\(u1-Tp*u0)

            x0 = unravel(x0, u11) 
            # compute+
        end


        if η<tol_η
            return (;solution=x0, message="Convergence", n_iterations=t)
        end

    end

    return (;message="No Convergence", n_iterations=T)

end



# function time_iteration(model;
#     T=500,
#     K=10,
#     tol_ε=1e-8,
#     tol_η=1e-6,
#     verbose=false
# )


#     N = length(model.grid)

#     x0 = GArray(model.grid, [SVector(model.x) for n=1:N])
#     x1 = deepcopy(x0)
#     dx = deepcopy(x0)

#     r0 = x0*0
    
#     J = dF0(model, x0, x0)[2]

#     local x0
#     local x1

#     for t=1:T
#         # r0 = F(model, x0, x0)
#         F!(r0, model, x0, x0)
#         ε = norm(r0)
#         if ε<tol_ε
#             break
#         end
#         if verbose
#             println("ϵ=$(ε)")
#         end
#         x1.data .= x0.data
#         for k=1:K
#             # r = F(model, x1, x0)
#             F!(r0, model, x1, x0)
#             # J = dF(model, x1, x0)
#             dF!(J, model, x1, x0)
#             # dx = J\r0
#             for n=1:length(r0)
#                 dx.data[n] = J.data[n]\r0.data[n]
#             end
#             e = norm(dx)
#             # println("e=$(e)")
#             x1.data .-= dx.data
#             if e<tol_η
#                 break
#             end
#         end
#         x0 = x1

#     end
#     return x0
# end



# # Alternative implementations

# function F0(model, s, x::SVector, xfut::GArray)
#     tot = SVector((x*0)...)
#     for (w, S) in τ(model, s, x)
#         ind = (S[1], S[3])
#         X = xfut(ind...)
#         tot += w*arbitrage(model,s,x,S,X)
#     end
#     return tot
# end


# F(model, controls::GArray, φ::GArray) =
#     GArray(
#         model.grid,
#         [F(model,s,x,φ) for (s,x) in zip(iti(model.grid), controls) ],
#     )

# function F!(out, model, controls, φ) 
#     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
#     n=0
#     for s in iti(model.grid)
#         n += 1
#         x = controls.data[n]
#         out.data[n] = F(model,s,x,φ)
#     end
#     # end
# end

# dF(model, controls::GArray, φ::GArray) =
#     GArray(    # this shouldn't be needed
#         model.grid,
#         [
#             ForwardDiff.jacobian(u->F(model, s, u, φ), x)
#             for (s,x) in zip(iti(model.grid), controls) 
#         ]
#     )

# function dF!(out, model, controls, φ) 
#     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
#     n=0
#     for s in iti(model.grid)
#         n += 1
#         x = controls.data[n]
#         out.data[n] = ForwardDiff.jacobian(u->F(model, s, u, φ), x)
#     end
#     # end
# end

# # function dF2!(out, model, controls, φ) 
# #     # for (n,(s,x)) in enumerate(zip(iti(model.grid), controls))
# #     n=0
# #     for s in iti(model.grid)
# #         n += 1
# #         x = controls.data[n]
# #         out.data[n] = ForwardDiff.jacobian(u->F(model, s, x, φ), φ)
# #     end
# #     # end
# # end


# FdF(model, controls::GArray, φ::GArray) =
#     GArray(
#         model.grid,
#         [
#             (F(model,s,x,φ), ForwardDiff.jacobian(u->F(model, s, u, φ), x))
#             for (s,x) in zip(iti(model.grid), controls) 
#         ]
#     )


# function F0(model, controls::GArray, xfut::GArray)

#     N = length(controls)
#     res = GArray(
#         model.grid,
#         zeros(typeof(controls[1]), N)
#     )
#     for (i,(s,x)) in enumerate(zip(iti(model.grid), controls))
#         res[i] = F(model,s,x,xfut)
#     end
#     return res
# end



# function dF0(model, controls::GArray, xfut::GArray)

#     N = length(controls)
#     res = deepcopy(controls)
#     dres = GArray(
#         model.grid,
#         zeros(typeof(res[1]*res[1]'), N)
#     )
#     for (i,(s,x)) in enumerate(zip(iti(model.grid), controls))
#         res[i] = F(model,s,x,xfut)
#         dres[i] = ForwardDiff.jacobian(u->F(model, s, u, xfut), x)
#     end
#     return res, dres
# end





