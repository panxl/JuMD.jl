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
    neighbor_list::NeighborList
    scache::Cache{Float64}
    vcache::Cache{SVector{D, Float64}}
end

function MMSystem(box, positions, masses, atomic_numbers, force_groups::ForceGroups{K,V}) where {K, V}
    N = Threads.nthreads()
    natoms = length(positions)
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
        cutoff = nothing
    elseif all(map(cutoff -> cutoff == first(cutoffs), cutoffs))
        cutoff = first(cutoffs)
        cell_list = LinkedCellList(natoms, cutoff, box, ratio=0.5)
        update!(cell_list, positions, box)
    else
        error("Cutoffs have to be the same for all force groups")
    end

    # check exclusion lists for force groups
    exclusions = []
    for group in force_groups.groups
        if hasfield(typeof(group), :exclusion) && !isnothing(group.exclusion)
            push!(exclusions, group.exclusion)
        end
    end
    if isempty(exclusions)
        exclusion = nothing
    elseif all(map(exclusion -> exclusion == first(exclusions), exclusions))
        exclusion = first(exclusions)
    else
        error("Exclusion lists have to be the same for all force groups")
    end

    # hard code maxnb and rskin for now
    maxnb = 1000
    rskin = 0.2
    neighbor_list = NeighborList(natoms, maxnb)
    if !isnothing(cutoff)
        rsch = cutoff + rskin
    else
        rsch = nothing
    end
    update!(neighbor_list, positions, box, rsch, exclusion, cell_list)

    D = length(eltype(positions))
    C = typeof(cell_list)

    MMSystem{D, N, C, K, V}(box, positions, velocities, forces, masses, atomic_numbers, force_groups, cell_list, neighbor_list, Cache(Float64), Cache(SVector{D, Float64}))
end

bounding_box(s::AbstractSystem)  = s.box
cell_list(s::AbstractSystem)     = s.cell_list
neighbor_list(s::AbstractSystem) = s.neighbor_list
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
