abstract type AbstractForce end

Base.@kwdef struct HarmonicBondForce <: AbstractForce
    indices::Vector{Tuple{Int, Int}} = []
    length::Vector{Float64} = []
    k::Vector{Float64} = []
end

function force!(system, f::HarmonicBondForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.length, f.k)
        i1, i2 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        r₀= f.length[i]
        k = f.k[i]

        (f1, f2), e = harmonic_bond_force(x1, x2, r₀, k)

        forces[i1] += f1
        forces[i2] += f2
        e_sum += e
    end
    return e_sum
end

function harmonic_bond_force(x1, x2, r₀, k)
    v = x2 - x1
    r = norm(v)
    Δr = r - r₀
    f1 = k * Δr / r * v
    f2 = -f1
    e = 0.5 * k * Δr^2
    return (f1, f2), e
end

Base.@kwdef struct HarmonicAngleForce <: AbstractForce
    indices::Vector{Tuple{Int, Int, Int}} = []
    angle::Vector{Float64} = []
    k::Vector{Float64} = []
end

function force!(system, f::HarmonicAngleForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.angle, f.k)
        i1, i2, i3 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        x3 = positions[i3]
        θ₀= f.angle[i]
        k = f.k[i]

        (f1, f2, f3), e = harmonic_angle_force(x1, x2, x3, θ₀, k)

        forces[i1] += f1
        forces[i2] += f2
        forces[i3] += f3
        e_sum += e
    end
    return e_sum
end

function harmonic_angle_force(x1, x2, x3, θ₀, k)
    v1 = x2 - x1
    v2 = x2 - x3
    r1 = norm(v1)
    r2 = norm(v2)
    cosine = clamp((v1 ⋅ v2) / (r1 * r2), -1, 1)
    θ = acos(cosine)
    Δθ = θ - θ₀
    dEdθ = k * Δθ

    cp = v1 × v2
    rp = norm(cp)
    f1 = dEdθ / (r1 * r1 * rp) * (v1 × cp)
    f3 = dEdθ / (r2 * r2 * rp) * (cp × v2)
    f2 = -(f1 + f3)
    e = 0.5 * k * Δθ^2
    return (f1, f2, f3), e
end

Base.@kwdef struct PeriodicTorsionForce <: AbstractForce
    indices::Vector{Tuple{Int, Int, Int, Int}} = []
    periodicity::Vector{Int} = []
    phase::Vector{Float64} = []
    k::Vector{Float64} = []
end

function force!(system, f::PeriodicTorsionForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.periodicity, f.phase, f.k)
        i1, i2, i3, i4 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        x3 = positions[i3]
        x4 = positions[i4]
        n = f.periodicity[i]
        ϕ₀ = f.phase[i]
        k = f.k[i]

        (f1, f2, f3, f4), e = periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)

        forces[i1] += f1
        forces[i2] += f2
        forces[i3] += f3
        forces[i4] += f4
        e_sum += e
    end
    return e_sum
end

function periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)
    v1 = x2 - x1
    v2 = x2 - x3
    v3 = x4 - x3

    r₂² = v2 ⋅ v2
    r₂ = sqrt(r₂²)
    cp1 = v1 × v2
    cp2 = v2 × v3
    x = cp1 ⋅ cp2
    y = (cp1 × cp2) ⋅ v2 / r₂
    ϕ = atan(y, x)

    Δϕ = n * ϕ - ϕ₀
    dEdϕ = -k * n * sin(Δϕ)

    f1 =  dEdϕ * r₂ / (cp1 ⋅ cp1) * cp1
    f4 = -dEdϕ * r₂ / (cp2 ⋅ cp2) * cp2

    s = (v1 ⋅ v2) / r₂² * f1 - (v3 ⋅ v2) / r₂² * f4
    f2 =  s - f1
    f3 = -s - f4

    e = k * (1.0 + cos(Δϕ))
    return (f1, f2, f3, f4), e
end

Base.@kwdef struct LennardJonesForce{C<:Union{Float64, Nothing}} <: AbstractForce
    sigma::Vector{Float64} = Float64[]
    epsilon::Vector{Float64} = Float64[]
    exclusion::Vector{Vector{Int}} = [Vector{Int}() for _ in 1 : length(sigma)]
    cutoff::C = nothing
end

function force!(system, forces::LennardJonesForce)
    e = force!(system, forces, system.cell_list)
    return e
end

function force!(system, f::LennardJonesForce, cl::NullCellList)
    natoms = length(f.epsilon)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    for i in 1 : (natoms - 1)
        # skip to next i-atom if i-atom's epsilon is zero
        if f.epsilon[i] == 0.0
            continue
        end

        for j in (i + 1 : natoms)
            # skip to next j-atom if j-atom's epsilon is zero
            if f.epsilon[j] == 0.0
                continue
            end

            # skip if j-atom is in i-atom's exclusion list
            if j ∈ f.exclusion[i]
                continue
            end

            x1 = positions[i]
            x2 = positions[j]
            v = x2 - x1

            σ = f.sigma[i] + f.sigma[j]
            ϵ = f.epsilon[i] * f.epsilon[j]

            e, ∂e∂v = lennard_jones_force(v, σ, ϵ)

            forces[i] += ∂e∂v
            forces[j] -= ∂e∂v
            e_sum += e
        end
    end
    return e_sum
end

function force!(system, f::LennardJonesForce, cl::LinkedCellList)
    positions = position(system)
    ncells = size(cl.head)
    rcut = f.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())

        @inbounds while i != 0
            # skip to next i-atom if i-atom's epsilon is zero
            if f.epsilon[i] == 0.0
                i = cl.list[i]
                continue
            end

            for offset in cl.offsets
                if offset == last(cl.offsets)
                    j = cl.list[i]
                else
                    cj = Tuple(ci + offset)
                    cj = mod.(cj, ncells)
                    j = cl.head[CartesianIndex(cj)]
                end

                while j != 0
                    # skip to next j-atom if j-atom's epsilon is zero
                    if f.epsilon[j] == 0.0
                        j = cl.list[j]
                        continue
                    end

                    # skip if j-atom is in i-atom's exclusion list
                    if j ∈ f.exclusion[i]
                        j = cl.list[j]
                        continue
                    end

                    x1 = positions[i]
                    x2 = positions[j]
                    v = x2 - x1

                    # apply minimum image convention
                    if any(system.box .!= 0.0)
                        v = minimum_image(v ./ system.box) .* system.box
                    end

                    # skip to next j-atom if r > rcut
                    r² = v ⋅ v
                    if r² > rcut²
                        j = cl.list[j]
                        continue
                    end

                    σ = f.sigma[i] + f.sigma[j]
                    ϵ = f.epsilon[i] * f.epsilon[j]

                    e, ∂e∂v = lennard_jones_force(v, σ, ϵ)

                    # apply switching function
                    r = sqrt(r²)
                    s, dsdr = shift(r, rcut)
                    ∂e∂v = ∂e∂v .* s + e * dsdr / r * v
                    e *= s

                    forces[i] += ∂e∂v
                    forces[j] -= ∂e∂v
                    e_threads[Threads.threadid()] += e

                    # Next j-atom
                    j = cl.list[j]
                end
            end

             # Next i-atom
            i = cl.list[i]
        end
    end # End loop over all cells

    e_sum = sum(e_threads)
    return e_sum
end

Base.@kwdef struct LennardJonesExceptionForce <: AbstractForce
    indices::Vector{Tuple{Int, Int}} = []
    sigma::Vector{Float64} = []
    epsilon::Vector{Float64} = []
end

function force!(system, f::LennardJonesExceptionForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.sigma, f.epsilon)
        i1, i2 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        v = x2 - x1

        # apply minimum image convention
        if any(system.box .!= 0.0)
            v = minimum_image(v ./ system.box) .* system.box
        end

        σ = f.sigma[i]
        ϵ = f.epsilon[i]

        e, ∂e∂v = lennard_jones_force(v, σ, ϵ)

        forces[i1] += ∂e∂v
        forces[i2] -= ∂e∂v
        e_sum += e
    end
    return e_sum
end

function lennard_jones_force(v, σ, ϵ)
    r² = v ⋅ v
    r = sqrt(r²)

    σ² = (σ / r)^2
    σ⁶ = σ² * σ² * σ²

    e = ϵ * σ⁶ * (σ⁶ - 1.0)
    ∂e∂v = -ϵ * σ⁶ * (12 * σ⁶ - 6.0) / r² * v
    return e, ∂e∂v
end

Base.@kwdef struct CoulombForce{C<:Union{Float64, Nothing}, E<:Union{<:AbstractRecip, Nothing}} <: AbstractForce
    charges::Vector{Float64} = Float64[]
    exclusion::Vector{Vector{Int}} = [Vector{Int}() for _ in 1 : length(charges)]
    cutoff::C = nothing
    recip::E = nothing
end

function force!(system, forces::CoulombForce)
    if isnothing(forces.recip)
        e = force!(system, forces, system.cell_list)
    else
        e = force!(system, forces, system.cell_list, forces.recip)
    end
    return e
end

function force!(system, f::CoulombForce, cl::NullCellList)
    natoms = length(f.charges)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    for i in 1 : (natoms - 1)
        # skip to next i-atom if i-atom's charge is zero
        if f.charges[i] == 0.0
            continue
        end

        for j in (i + 1) : natoms
            # skip to next j-atom if j-atom's charge is zero
            if f.charges[j] == 0.0
                continue
            end

            # skip if j-atom is in i-atom's exclusion list
            if j ∈ f.exclusion[i]
                continue
            end

            x1 = positions[i]
            x2 = positions[j]
            v = x2 - x1

            fac = KE * f.charges[i] * f.charges[j]
            e, ∂e∂v = fac .* coulomb_potential(v)

            forces[i] += ∂e∂v
            forces[j] -= ∂e∂v
            e_sum += e
        end
    end
    return e_sum
end

function force!(system, f::CoulombForce, cl::LinkedCellList)
    positions = position(system)
    ncells = size(cl.head)
    rcut = f.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())

        @inbounds while i != 0
            # skip to next i-atom if i-atom's charge is zero
            if f.charges[i] == 0.0
                i = cl.list[i]
                continue
            end

            for offset in cl.offsets
                if offset == last(cl.offsets)
                    j = cl.list[i]
                else
                    cj = Tuple(ci + offset)
                    cj = mod.(cj, ncells)
                    j = cl.head[CartesianIndex(cj)]
                end

                while j != 0
                    # skip to next j-atom if j-atom's charge is zero
                    if f.charges[j] == 0.0
                        j = cl.list[j]
                        continue
                    end

                    # skip if j-atom is in i-atom's exclusion list
                    if j ∈ f.exclusion[i]
                        j = cl.list[j]
                        continue
                    end

                    x1 = positions[i]
                    x2 = positions[j]
                    v = x2 - x1

                    # apply minimum image convention
                    if any(system.box .!= 0.0)
                        v = minimum_image(v ./ system.box) .* system.box
                    end

                    # skip to next j-atom if r > rcut
                    r² = v ⋅ v
                    if r² > rcut²
                        j = cl.list[j]
                        continue
                    end

                    fac = KE * f.charges[i] * f.charges[j]
                    e, ∂e∂v = fac .* coulomb_potential(v)

                    # apply switching function
                    r = sqrt(r²)
                    s, dsdr = shift(r, rcut)
                    ∂e∂v = ∂e∂v .* s + e * dsdr / r * v
                    e *= s

                    forces[i] += ∂e∂v
                    forces[j] -= ∂e∂v
                    e_threads[Threads.threadid()] += e

                    # Next j-atom
                    j = cl.list[j]
                end
            end

            # Next i-atom
            i = cl.list[i]
        end
    end # End loop over all cells

    e_sum = sum(e_threads)
    return e_sum
end

Base.@kwdef struct CoulombExceptionForce <: AbstractForce
    indices::Vector{Tuple{Int, Int}} = []
    charge_prod::Vector{Float64} = []
end

function force!(system, f::CoulombExceptionForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.charge_prod)
        i1, i2 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        v = x2 - x1

        # apply minimum image convention
        if any(system.box .!= 0.0)
            v = minimum_image(v ./ system.box) .* system.box
        end

        fac = KE * f.charge_prod[i]
        e, ∂e∂v = fac .* coulomb_potential(v)

        forces[i1] += ∂e∂v
        forces[i2] -= ∂e∂v
        e_sum += e
    end
    return e_sum
end

function coulomb_potential(v)
    r² = v ⋅ v
    r = sqrt(r²)
    e = 1 / r
    ∂e∂v = -e / r² * v
    return e, ∂e∂v
end

function force!(system, f::CoulombForce, cl::LinkedCellList, recip::AbstractRecip)
    positions = position(system)
    ncells = size(cl.head)
    rcut = f.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())

        @inbounds while i != 0
            # skip to next i-atom if i-atom's charge is zero
            if f.charges[i] == 0.0
                i = cl.list[i]
                continue
            end

            for offset in cl.offsets
                if offset == last(cl.offsets)
                    j = cl.list[i]
                else
                    cj = Tuple(ci + offset)
                    cj = mod.(cj, ncells)
                    j = cl.head[CartesianIndex(cj)]
                end

                while j != 0
                    # skip to next j-atom if j-atom's charge is zero
                    if f.charges[j] == 0.0
                        j = cl.list[j]
                        continue
                    end

                    x1 = positions[i]
                    x2 = positions[j]
                    v = x2 - x1

                    # apply minimum image convention
                    v = minimum_image(v ./ system.box) .* system.box

                    # depends on if j-atom is in i-atom's exclusion list,
                    # we either substract k-space or add real space
                    # contributions from the unit cell
                    if j ∈ f.exclusion[i]
                        fac = -KE * f.charges[i] * f.charges[j]
                        e, ∂e∂v = fac .* ewald_recip_potential(v, recip.alpha)
                    else
                        # skip to next j-atom if r > rcut
                        r² = v ⋅ v
                        if r² > rcut²
                            j = cl.list[j]
                            continue
                        end

                        fac = KE * f.charges[i] * f.charges[j]
                        e, ∂e∂v = fac .* ewald_real_potential(v, recip.alpha)
                    end

                    forces[i] += ∂e∂v
                    forces[j] -= ∂e∂v
                    e_threads[Threads.threadid()] += e

                    # Next j-atom
                    j = cl.list[j]
                end
            end

            # Next i-atom
            i = cl.list[i]
        end
    end # End loop over all cells

    e_sum = sum(e_threads)

    # k-space

    e_recip = recip(system.box, f.charges, positions, system.forces)

    e_sum += e_recip

    return e_sum
end

function ewald_real_potential(v, α)
    r² = v ⋅ v
    r = sqrt(r²)
    e = erfc(α * r) / r
    ∂e∂v = -(e + 2 * α * exp(-(α * r)^2) / SQRTPI) / r² * v
    return e, ∂e∂v
end

function ewald_recip_potential(v, α)
    r² = v ⋅ v
    r = sqrt(r²)
    e = erf(α * r) / r
    ∂e∂v = -(e - 2 * α * exp(-(α * r)^2) / SQRTPI) / r² * v
    return e, ∂e∂v
end
