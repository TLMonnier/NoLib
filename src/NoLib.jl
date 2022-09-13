module NoLib

    import Base: eltype
    using LabelledArrays
    using StaticArrays
    using ForwardDiff
    
    include("interp.jl")
    using NoLib.Interpolation: interp
    
    ⟂(a,b) = min(a,b)
    # ⟂ᶠ(a,b)

    import Base: getindex
    
    include("grids.jl")
    include("garray.jl")
    include("model.jl")
    include("simul.jl")
    include("time_iteration.jl")
    include("hetag.jl")

end # module

