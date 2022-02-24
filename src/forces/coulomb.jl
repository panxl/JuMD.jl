Base.@kwdef struct CoulombForce{C<:Union{Float64, Nothing}, R<:AbstractRecip} <: AbstractForce
    charges::Vector{Float64} = Float64[]
    exclusion::Vector{Vector{Int}} = [Vector{Int}() for _ in 1 : length(charges)]
    cutoff::C = nothing
    recip::R = NullRecip()
end

function force!(system, forces::CoulombForce)
    if cell_list(system) isa NullCellList
        e = force!(system, forces, cell_list(system), forces.recip)
    else
        e = force!(system, forces, neighbor_list(system), forces.recip)
    end
    return e
end

function force!(system, f::CoulombForce, cl::NullCellList, recip::NullRecip)
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

function force!(system, f::CoulombForce, cl::LinkedCellList, recip::NullRecip)
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

function force!(system, f::CoulombForce, nbl::NeighborList, recip::AbstractRecip)
    positions = position(system)
    rcut = f.cutoff
    rcut² = rcut^2
    e_threads = zeros(Threads.nthreads())

    Threads.@threads for i in 1 : length(positions)
        forces = force(system, Threads.threadid())

        # skip to next i-atom if i-atom's charge is zero
        q₁ = f.charges[i]
        if q₁ == 0.0
            continue
        end

        # load other i-atom's information
        x1 = positions[i]

        ilist = nbl.list[i]
        @inbounds for j_idx in 1 : nbl.n[i]
            j = ilist[j_idx]

            # skip to next j-atom if j-atom's epsilon is zero
            q₂ = f.charges[j]
            if q₂ == 0.0
                continue
            end

            # load other j-atom's information
            x2 = positions[j]

            v = x2 - x1

            # apply minimum image convention
            v = minimum_image(v ./ system.box) .* system.box

            # skip to next j-atom if r > rcut
            r² = v ⋅ v
            if r² > rcut²
                continue
            end

            fac = KE * q₁ * q₂
            e, ∂e∂v = fac .* ewald_real_potential(v, recip.alpha)

            forces[i] += ∂e∂v
            forces[j] -= ∂e∂v
            e_threads[Threads.threadid()] += e
        end

        #substract k-space contributions from the unit cell
        @inbounds for j in f.exclusion[i]
            if j > i
                continue
            end

            x2 = positions[j]
            v = x2 - x1

            # apply minimum image convention
            v = minimum_image(v ./ system.box) .* system.box

            q₂ = f.charges[j]

            fac = -KE * q₁ * q₂
            e, ∂e∂v = fac .* ewald_recip_potential(v, recip.alpha)

            forces[i] += ∂e∂v
            forces[j] -= ∂e∂v
            e_threads[Threads.threadid()] += e
        end
    end

    e_sum = sum(e_threads)

    # k-space

    e_recip = recip(system.box, f.charges, positions, system.forces)

    e_sum += e_recip

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
