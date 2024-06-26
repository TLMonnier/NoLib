
### transition function


function τ(model, ss::Tuple, a::SVector)


    p = model.p
    (i,_),(m, s) = ss # get current state values

    Q = model.grid.g1.points

    j = 1

    it = (
        (
            model.P[i,j],
            (
                (j,),
                (
                    Q[j],
                    transition(model, m, s, a, Q[j], p)
                )
            )
        )
        for j in 1:size(model.P, 2)
    )

    it

end


function trembling__hand(g::CGrid{1}, xv)
    x = xv[1]
    r = g.ranges[1]
    u = (x-r[1])/(r[2]-r[1])
    n = r[3]

    i_ = floor(Int, u*(n-1))
    i_ = max(min(i_,n-2),0)
    λ = u*(n-1)-i_

    λ = min(max(λ, 0.0), 1.0)
    
    (
        (1-λ, i_+1),
        (λ, i_+2)
    )
end

using ResumableFunctions

@resumable function τ_fit(model, ss::Tuple, a::SVector; linear_index=false)

    p = model.p

    i = ss[1][1]

    Q = model.grid.g1.points

    n_m = length(model.m)

    # (m, s) = (ss[2][1:n_m], ss[2][n_m+1:end])
    # (m, s) = (ss[2][1], ss[2][2])
    (i,_),(m, s) = ss


    for j in 1:size(model.P, 2)

        S = transition(model, m, s, a, Q[j], p)

        for (w, i_S) in trembling__hand(model.grid.g2, S)

            res = (
                model.P[i,j]*w,

                (
                    (linear_index ? to__linear_index(model.grid, (j,i_S)) : (j,i_S)),

                    (Q[j], model.grid.g2[i_S])
                )
            )
            if linear_index
                res
                # res::Tuple{Float64, Tuple{Int64, Tuple{SVector{1},SVector{1}}}}
            else
                res
                # res::Tuple{Float64, Tuple{Tuple{Int64, Int64}, Tuple{SVector{1},SVector{1}}}}
            end
            @yield res
        end
    end
end

function G(model, μ::GDist{T}, x) where T
    μ1 = GArray(μ.grid, zeros(Float64, length(μ)))
    for ss in iti(model.grid)
        a = x[ss[1]...]
        for (w, (ind, _)) in τ_fit(model, ss, a)
            μ1[ind...] += w*μ[ind...]
        end
    end

    μ1

end

# TODO: write interpolation version of G


function transition_matrix(model, x)

    N = length(x)
    P = zeros(N,N)

    for (ss,a) in zip(iti(model.grid),x)
        ind_i = ss[1]
        i = to__linear_index(model.grid, ind_i)
        for (w, (ind_j, _)) in τ_fit(model, ss, a)
            j = to__linear_index(model.grid, ind_j)
            P[i,j] = w
        end
    end

    P

end


using LinearAlgebra: I

function ergodic_distribution(model, x::GArray)

    P = transition_matrix(model, x)

    PP = [ (P-I)[:,1:end-1] ;;  ones(size(P,1))]
    R = zeros(size(PP,1))
    R[end] = 1
    μ = PP'\R
    
    ergo = GArray(model.grid, μ)

    return ergo

end


## TODO: some ways to plot the ergo dist...

τ(model, ss::Tuple, φ) = τ(model, ss, φ(ss))

## TODO: some simulation

# TODO
# function simulate(model, ss::Tuple, φ; T=20)
#     Typ = typeof(ss)
#     sim = Typ[ss]
#     for t = 1:T
#     end
# end