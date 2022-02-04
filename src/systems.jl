struct MMSystem{D, K, V} <:AbstractSystem{D}
    positions::Vector{SVector{D, Float64}}
    velocities::Vector{SVector{D, Float64}}
    forces::Vector{SVector{D, Float64}}
    masses::Vector{Float64}
    atomic_numbers::Vector{Int}
    force_groups::ForceGroups{K, V}
    scache::Cache{Float64}
    vcache::Cache{SVector{D, Float64}}
end

function MMSystem(positions::AbstractVector{SVector{D, F}},
                velocities::AbstractVector{SVector{D, F}},
                forces::AbstractVector{SVector{D, F}},
                masses::AbstractVector{F},
                atomic_numbers::AbstractVector{I},
                force_groups::ForceGroups{K,V}) where {D, F<:AbstractFloat, I<:Integer, K, V}
    MMSystem{D, K, V}(positions, velocities, forces, masses, atomic_numbers, force_groups, Cache(Float64), Cache(SVector{D, Float64}))
end

position(s::AbstractSystem)      = s.positions
velocity(s::AbstractSystem)      = s.velocities
force(s::AbstractSystem)         = s.forces
atomic_number(s::AbstractSystem) = s.atomic_numbers
mass(s::AbstractSystem)          = s.masses

inverse_mass(s::AbstractSystem) = !hasproperty(s.scache, :inverse_masses) ? s.scache.inverse_masses = 1 ./ s.masses : s.scache.inverse_masses
inverse_mass!(s::AbstractSystem) = hasproperty(s.scache, :inverse_masses) ? s.scache.inverse_masses .= 1 ./ s.masses : s.scache.inverse_masses = 1 ./ s.masses

velocity_half(s::AbstractSystem) = !hasproperty(s.vcache, :velocities_half) ? s.vcache.velocities_half = similar(s.velocities) : s.vcache.velocities_half
position_last(s::AbstractSystem) = !hasproperty(s.vcache, :positions_last) ? s.vcache.positions_last = similar(s.positions) : s.vcache.positions_last

function force!(s::AbstractSystem, force_group::ForceGroup{<:AbstractForce})
    e = 0.0
    for f in force_group.forces
        e += force!(s, f)
    end
    return e
end

function force!(s::AbstractSystem, force_groups::ForceGroups)
    fill!(s.forces, zeros(eltype(s.forces)))
    force_groups.energies .= map(fg::ForceGroup -> force!(s, fg), values(force_groups.groups))
    return nothing
end

force!(s::AbstractSystem) = force!(s, s.force_groups)