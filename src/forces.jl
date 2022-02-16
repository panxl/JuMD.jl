abstract type AbstractSystem{D} end

abstract type AbstractForce end

struct HarmonicBondForce <: AbstractForce
    indices::Tuple{Int, Int}
    length::Float64
    k::Float64
end

function force!(system, f::HarmonicBondForce)
    positions = position(system)
    forces = force(system)
    i1, i2 = f.indices

    x1 = positions[i1]
    x2 = positions[i2]
    r₀= f.length
    k = f.k

    (f1, f2), e = harmonic_bond_force(x1, x2, r₀, k)
    forces[i1] += f1
    forces[i2] += f2

    return e
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

struct HarmonicAngleForce <: AbstractForce
    indices::Tuple{Int, Int, Int}
    angle::Float64
    k::Float64
end

function force!(system, f::HarmonicAngleForce)
    positions = position(system)
    forces = force(system)
    i1, i2, i3 = f.indices

    x1 = positions[i1]
    x2 = positions[i2]
    x3 = positions[i3]
    θ₀= f.angle
    k = f.k

    (f1, f2, f3), e = harmonic_angle_force(x1, x2, x3, θ₀, k)
    forces[i1] += f1
    forces[i2] += f2
    forces[i3] += f3

    return e
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

struct PeriodicTorsionForce <: AbstractForce
    indices::Tuple{Int, Int, Int, Int}
    periodicity::Int
    phase::Float64
    k::Float64
end

function force!(system, f::PeriodicTorsionForce)
    positions = position(system)
    forces = force(system)
    i1, i2, i3, i4 = f.indices

    x1 = positions[i1]
    x2 = positions[i2]
    x3 = positions[i3]
    x4 = positions[i4]
    n = f.periodicity
    ϕ₀ = f.phase
    k = f.k

    (f1, f2, f3, f4), e = periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)
    forces[i1] += f1
    forces[i2] += f2
    forces[i3] += f3
    forces[i4] += f4

    return e
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

struct LennardJonesForce <: AbstractForce
    sigma::Vector{Float64}
    epsilon::Vector{Float64}
    exlusion::Vector{Set{Int}}
end

LennardJonesForce(sigma, epsilon) = LennardJonesForce(sigma, epsilon, [Set{Int}() for _ in 1 : length(sigma)])

function force!(system::AbstractSystem, f::LennardJonesForce, cl::NullCellList)
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

            # skip if j-atom is in i-atom's exlusion list
            if j ∈ f.exlusion[i]
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

function force!(system::AbstractSystem, f::LennardJonesForce, cl::LinkedCellList)
    positions = position(system)
    ncells = size(cl.head)
    rcut = system.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())

        while i != 0
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

                    # skip if j-atom is in i-atom's exlusion list
                    if j ∈ f.exlusion[i]
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

struct LennardJonesExceptionForce <: AbstractForce
    indices::Tuple{Int, Int}
    sigma::Float64
    epsilon::Float64
end

function force!(system::AbstractSystem, f::LennardJonesExceptionForce)
    positions = position(system)
    forces = force(system)
    i, j = f.indices

    x1 = positions[i]
    x2 = positions[j]
    v = x2 - x1

    # apply minimum image convention
    if any(system.box .!= 0.0)
        v = minimum_image(v ./ system.box) .* system.box
    end

    σ = f.sigma
    ϵ = f.epsilon

    e, ∂e∂v = lennard_jones_force(v, σ, ϵ)
    forces[i] += ∂e∂v
    forces[j] -= ∂e∂v

    return e
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

struct CoulombForce{E<:Union{Ewald,Nothing}} <: AbstractForce
    charge::Vector{Float64}
    exlusion::Vector{Set{Int}}
    ewald::E
end

CoulombForce(charge) = CoulombForce(charge, [Set{Int}() for _ in 1 : length(charge)], nothing)

function force!(system::AbstractSystem, f::CoulombForce, cl::NullCellList)
    natoms = length(f.charge)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    for i in 1 : (natoms - 1)
        # skip to next i-atom if i-atom's charge is zero
        if f.charge[i] == 0.0
            continue
        end

        for j in (i + 1) : natoms
            # skip to next j-atom if j-atom's charge is zero
            if f.charge[j] == 0.0
                continue
            end

            # skip if j-atom is in i-atom's exlusion list
            if j ∈ f.exlusion[i]
                continue
            end

            x1 = positions[i]
            x2 = positions[j]
            v = x2 - x1

            fac = KE * f.charge[i] * f.charge[j]
            e, ∂e∂v = fac .* coulomb_potential(v)

            forces[i] += ∂e∂v
            forces[j] -= ∂e∂v
            e_sum += e
        end
    end
    return e_sum
end

function force!(system::AbstractSystem, f::CoulombForce, cl::LinkedCellList)
    positions = position(system)
    ncells = size(cl.head)
    rcut = system.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())

        while i != 0
            # skip to next i-atom if i-atom's charge is zero
            if f.charge[i] == 0.0
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
                    if f.charge[j] == 0.0
                        j = cl.list[j]
                        continue
                    end

                    # skip if j-atom is in i-atom's exlusion list
                    if j ∈ f.exlusion[i]
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

                    fac = KE * f.charge[i] * f.charge[j]
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

struct CoulombExceptionForce <: AbstractForce
    indices::Tuple{Int, Int}
    charges::Tuple{Float64, Float64}
    scaling::Float64
end

function force!(system::AbstractSystem, f::CoulombExceptionForce)
    positions = position(system)
    forces = force(system)
    i, j = f.indices

    x1 = positions[i]
    x2 = positions[j]
    v = x2 - x1

    # apply minimum image convention
    if any(system.box .!= 0.0)
        v = minimum_image(v ./ system.box) .* system.box
    end

    fac = f.scaling * KE * f.charges[1] * f.charges[2]
    e, ∂e∂v = fac .* coulomb_potential(v)

    forces[i] += ∂e∂v
    forces[j] -= ∂e∂v

    return e
end

function coulomb_potential(v)
    r² = v ⋅ v
    r = sqrt(r²)
    e = 1 / r
    ∂e∂v = -e / r² * v
    return e, ∂e∂v
end

function force!(system::AbstractSystem, f::CoulombForce, cl::LinkedCellList, ewald::Ewald)
    positions = position(system)
    ncells = size(cl.head)
    rcut = system.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())

        while i != 0
            # skip to next i-atom if i-atom's charge is zero
            if f.charge[i] == 0.0
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
                    if f.charge[j] == 0.0
                        j = cl.list[j]
                        continue
                    end

                    x1 = positions[i]
                    x2 = positions[j]
                    v = x2 - x1

                    # apply minimum image convention
                    v = minimum_image(v ./ system.box) .* system.box

                    # depends on if j-atom is in i-atom's exlusion list,
                    # we either substract k-space or add real space
                    # contributions from the unit cell
                    if j ∈ f.exlusion[i]
                        fac = -KE * f.charge[i] * f.charge[j]
                        e, ∂e∂v = fac .* ewald_reciporical_potential(v, ewald.alpha)
                    else
                        # skip to next j-atom if r > rcut
                        r² = v ⋅ v
                        if r² > rcut²
                            j = cl.list[j]
                            continue
                        end

                        fac = KE * f.charge[i] * f.charge[j]
                        e, ∂e∂v = fac .* ewald_real_potential(v, ewald.alpha)
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

    fractional_coordinates = fractional_coordinate(system)
    update!(ewald, fractional_coordinates)

    charges = f.charge
    recip_box = 2 * pi ./ system.box
    inv_v = 1 / prod(system.box)

    fill!(e_threads, zero(eltype(e_threads)))

    Threads.@threads for i in eachindex(ewald.kvectors, ewald.kfactors)
        forces = force(system, Threads.threadid())
        qeir = ewald._qeir[Threads.threadid()]
        kx, ky, kz = ewald.kvectors[i]
        kfac = ewald.kfactors[i]
        kr = ewald.kvectors[i] .* recip_box

        term = zero(ComplexF64)
        @inbounds @simd for j in 1 : length(charges)
            qeir[j] = charges[j] * ewald.eir[j, kx, 1] * ewald.eir[j, ky, 2] * ewald.eir[j, kz, 3]
            term += qeir[j]
        end
        e_threads[Threads.threadid()] += kfac * real(conj(term) * term)

        @inbounds @simd for j in 1 : length(charges)
            f_recip = 2 * kfac * imag(conj(term) * qeir[j]) .* kr
            forces[j] += f_recip .* (KE * inv_v)
        end
    end

    e_recip = sum(e_threads) * KE / prod(system.box)

    # substract self part of k-space sum
    e_recip -= KE * ewald.alpha * (charges ⋅ charges) / SQRTPI

    e_sum += e_recip

    return e_sum
end
