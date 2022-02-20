abstract type AbstractSystem end

struct MMSystem{D, N, C, K, V} <: AbstractSystem
    box::SVector{D, Float64}
    positions::Vector{SVector{D, Float64}}
    velocities::Vector{SVector{D, Float64}}
    forces::NTuple{N, Vector{SVector{D, Float64}}}
    masses::Vector{Float64}
    atomic_numbers::Vector{Int}
    force_groups::ForceGroups{K, V}
    cell_list::C
    scache::Cache{Float64}
    vcache::Cache{SVector{D, Float64}}
end

function MMSystem(box, positions, masses, atomic_numbers, force_groups::ForceGroups{K,V}) where {K, V}
    N = Threads.nthreads()
    positions = SVector.(positions)
    velocities = zero(positions)
    forces = ntuple(i -> zero(positions), N)

    # check cutoff for force groups and set cell list accordingly
    cutoffs = []
    for group in force_groups.groups
        if hasfield(typeof(group), :cutoff) && !isnothing(group.cutoff)
            push!(cutoffs, group.cutoff)
        end
    end
    if isempty(cutoffs)
        cell_list = NullCellList()
    elseif all(map(cutoff -> cutoff == first(cutoffs), cutoffs))
        cell_list = LinkedCellList(length(positions), first(cutoffs), box)
    else
        error("Cutoffs have to be the same for all force groups")
    end

    D = length(eltype(positions))
    C = typeof(cell_list)

    MMSystem{D, N, C, K, V}(box, positions, velocities, forces, masses, atomic_numbers, force_groups, cell_list, Cache(Float64), Cache(SVector{D, Float64}))
end

bounding_box(s::AbstractSystem)  = s.box
cell_list(s::AbstractSystem)     = s.cell_list
position(s::AbstractSystem)      = s.positions
velocity(s::AbstractSystem)      = s.velocities
force(s::AbstractSystem)         = s.forces[1]
force(s::AbstractSystem, i::Integer) = s.forces[i]
atomic_number(s::AbstractSystem) = s.atomic_numbers
mass(s::AbstractSystem)          = s.masses

inverse_mass(s::AbstractSystem)  = !hasproperty(s.scache, :inverse_masses) ? s.scache.inverse_masses = 1 ./ s.masses : s.scache.inverse_masses
inverse_mass!(s::AbstractSystem) = hasproperty(s.scache, :inverse_masses) ? s.scache.inverse_masses .= 1 ./ s.masses : s.scache.inverse_masses = 1 ./ s.masses

velocity_half(s::AbstractSystem) = !hasproperty(s.vcache, :velocities_half) ? s.vcache.velocities_half = similar(s.velocities) : s.vcache.velocities_half
position_last(s::AbstractSystem) = !hasproperty(s.vcache, :positions_last) ? s.vcache.positions_last = similar(s.positions) : s.vcache.positions_last
