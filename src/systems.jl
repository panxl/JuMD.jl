struct MMSystem{D, C, N, K, V} <:AbstractSystem{D}
    box::SVector{D, Float64}
    positions::Vector{SVector{D, Float64}}
    velocities::Vector{SVector{D, Float64}}
    forces::NTuple{N, Vector{SVector{D, Float64}}}
    masses::Vector{Float64}
    atomic_numbers::Vector{Int}
    force_groups::ForceGroups{K, V}
    cell_list::C
    cutoff::Float64
    scache::Cache{Float64}
    vcache::Cache{SVector{D, Float64}}
end

function MMSystem(box,
                  positions,
                  masses,
                  atomic_numbers,
                  force_groups::ForceGroups{K,V};
                  cutoff=nothing) where {K, V}
    positions = SVector.(positions)
    N = Threads.nthreads()
    velocities = zero(positions)
    forces = ntuple(i -> zero(positions), N)

    if isnothing(cutoff)
        cell_list = NullCellList()
        rcut = 0.0
    else
        rcut = float(cutoff)
        cell_list = LinkedCellList(length(positions), rcut, box)
    end
    D = length(eltype(positions))
    C = typeof(cell_list)

    MMSystem{D, C, N, K, V}(box, positions, velocities, forces, masses, atomic_numbers, force_groups, cell_list, rcut, Cache(Float64), Cache(SVector{D, Float64}))
end

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

function fractional_coordinate(s::AbstractSystem)
    if !hasproperty(s.vcache, :fractional_coordinates)
        s.vcache.fractional_coordinates = similar(s.positions)
        for i in eachindex(s.positions)
            s.vcache.fractional_coordinates[i] = s.positions[i] ./ s.box
        end
    end
    return s.vcache.fractional_coordinates
end

function fractional_coordinate!(s::AbstractSystem)
    if !hasproperty(s.vcache, :fractional_coordinates)
        s.vcache.fractional_coordinates = similar(s.positions)
    end
    for i in eachindex(s.positions)
        s.vcache.fractional_coordinates[i] = s.positions[i] ./ s.box
    end
    return s.vcache.fractional_coordinates
end

function force!(s::AbstractSystem, force_group::ForceGroup{<:Vector{<:AbstractForce}})
    e = 0.0
    for f in force_group.forces
        e += force!(s, f)
    end
    return e
end

function force!(s::AbstractSystem, force_group::ForceGroup{<:CoulombForce})
    if isnothing(force_group.forces.recip)
        e = force!(s, force_group.forces, s.cell_list)
    else
        e = force!(s, force_group.forces, s.cell_list, force_group.forces.recip)
    end
    return e
end

function force!(s::AbstractSystem, force_group::ForceGroup{<:AbstractForce})
    e = force!(s, force_group.forces, s.cell_list)
    return e
end

function force!(s::AbstractSystem, force_groups::ForceGroups)
    for i in 1:Threads.nthreads()
        forces = force(s, i)
        fill!(forces, zeros(eltype(forces)))
    end
    r_box = fractional_coordinate!(s)
    update!(s.cell_list, r_box)
    force_groups.energies .= map(fg::ForceGroup -> force!(s, fg), values(force_groups.groups))
    for i in 2:Threads.nthreads()
        forces = force(s, i)
        force(s) .+= forces
        fill!(forces, zeros(eltype(forces)))
    end
    return nothing
end

force!(s::AbstractSystem) = force!(s, s.force_groups)
