abstract type AbstractForce end

include("forces/bonded.jl")
include("forces/lennard_jones.jl")
include("forces/ewald.jl")
include("forces/pme.jl")
include("forces/coulomb.jl")

Base.Base.@kwdef struct ForceGroups{K, V}
    groups::NamedTuple{K, V}
    energies::Vector{Float64} = zeros(Float64, length(groups))
    nrespa::Vector{Int} = ones(Int, length(groups))
end

function force!(system, force_groups::ForceGroups)
    for i in 1:Threads.nthreads()
        forces = force(system, i)
        fill!(forces, zeros(eltype(forces)))
    end
    update!(cell_list(system), position(system), bounding_box(system))
    force_groups.energies .= map(f -> force!(system, f), values(force_groups.groups))
    for i in 2:Threads.nthreads()
        forces = force(system, i)
        force(system) .+= forces
        fill!(forces, zeros(eltype(forces)))
    end
    return nothing
end

force!(system) = force!(system, system.force_groups)
