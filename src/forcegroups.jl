struct ForceGroup{F<:AbstractForce, N<:Union{Float64,Vector{Float64}}}
    forces::Vector{F}
    scaling_factors::N
end

ForceGroup(forces) = ForceGroup(forces, 1.0)
ForceGroup(F::Type) = ForceGroup(Vector{F}())

struct ForceGroups{K, V}
    groups::NamedTuple{K, V}
    energies::Vector{Float64}
    nrespa::Vector{Int}
end

ForceGroups(groups) = ForceGroups(groups, zeros(Float64, length(groups)), ones(Int, length(groups)))