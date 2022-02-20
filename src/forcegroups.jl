Base.Base.@kwdef struct ForceGroups{K, V}
    groups::NamedTuple{K, V}
    energies::Vector{Float64} = zeros(Float64, length(groups))
    nrespa::Vector{Int} = ones(Int, length(groups))
end
