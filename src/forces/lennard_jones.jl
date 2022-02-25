Base.@kwdef struct LennardJonesForce{C<:Union{Float64, Nothing}} <: AbstractForce
    epsilon::Vector{Float64} = Float64[]
    sigma::Vector{Float64} = Float64[]
    exclusion::Vector{Vector{Int}} = [Vector{Int}() for _ in 1 : length(epsilon)]
    cutoff::C = nothing
end

function force!(system, forces::LennardJonesForce)
    if cell_list(system) isa NullCellList
        e = force!(system, forces, cell_list(system))
    else
        e = force!(system, forces, neighbor_list(system), system.soa)
    end
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

        for j in (i + 1) : natoms
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

    @batch for ci in CartesianIndices(cl.head)
        i = cl.head[ci]
        forces = force(system, Threads.threadid())
        e_thread = 0.0

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
                    e_thread += e

                    # Next j-atom
                    j = cl.list[j]
                end
            end

            # Next i-atom
            i = cl.list[i]
        end

        e_threads[Threads.threadid()] += e_thread
    end # End loop over all cells

    e_sum = sum(e_threads)
    return e_sum
end

function force!(system, f::LennardJonesForce, nbl::NeighborList)
    positions = position(system)
    rcut = f.cutoff
    rcut² = rcut^2
    @unpack n_neighbors, lists = nbl

    e_threads = zeros(Threads.nthreads())

    @batch for i in 1 : length(positions)
        forces = force(system, Threads.threadid())

        # skip to next i-atom if i-atom's epsilon is zero
        ϵ₁ = f.epsilon[i]
        if ϵ₁ == 0.0
            continue
        end

        # load i-atom's other information
        x1 = positions[i]
        σ₁ = f.sigma[i]

        # accumulator for i-atom
        e_thread = 0.0
        f_thread = zeros(SVector{3, Float64})

        list = lists[i]
        n = n_neighbors[i]

        for j_idx in 1 : n
            j = list[j_idx]
            ϵ₂ = f.epsilon[j]

            # skip to next j-atom if j-atom's epsilon is zero
            if ϵ₂ == 0.0
                continue
            end

            # load j-atom's other information
            x2 = positions[j]
            σ₂ = f.sigma[j]

            v = x2 - x1

            # apply minimum image convention
            if any(system.box .!= 0.0)
                v = minimum_image(v ./ system.box) .* system.box
            end

            # skip to next j-atom if r > rcut
            r² = v ⋅ v
            if r² > rcut²
                continue
            end

            σ = σ₁ + σ₂
            ϵ = ϵ₁ * ϵ₂

            e, ∂e∂v = lennard_jones_force(v, σ, ϵ)

            # apply switching function
            r = sqrt(r²)
            s, dsdr = shift(r, rcut)
            ∂e∂v = ∂e∂v .* s + e * dsdr / r * v
            e *= s

            # accumulate
            e_thread += e
            f_thread += ∂e∂v
            forces[j] -= ∂e∂v
        end

        forces[i] += f_thread
        e_threads[Threads.threadid()] += e_thread
    end

    e_sum = sum(e_threads)
    return e_sum
end

function force!(system, f::LennardJonesForce, nbl::NeighborList, soa)
    rcut = f.cutoff
    rcut² = rcut^2

    @unpack n_neighbors, lists = nbl
    @unpack x, y, z, fx, fy, fz = soa

    e_threads = zeros(Threads.nthreads())

    for i in 1 : length(system.positions)
        ϵ₁ = f.epsilon[i]

        # skip to next i-atom if i-atom's epsilon is zero
        if ϵ₁ == 0.0
            continue
        end

        # load i-atom's other information
        x1 = x[i]
        y1 = y[i]
        z1 = z[i]
        σ₁ = f.sigma[i]

        # accumulator for i-atom
        e_thread = fx_thread = fy_thread = fz_thread = 0.0

        list = lists[i]
        n = n_neighbors[i]

        @turbo for j_idx in 1 : n
            j = list[j_idx]

            # load j-atom's information
            x2 = x[j]
            y2 = y[j]
            z2 = z[j]
            ϵ₂ = f.epsilon[j]
            σ₂ = f.sigma[j]

            vx = x2 - x1
            vy = y2 - y1
            vz = z2 - z1

            # apply minimum image convention
            vx = minimum_image(vx / system.box[1]) * system.box[1]
            vy = minimum_image(vy / system.box[2]) * system.box[2]
            vz = minimum_image(vz / system.box[3]) * system.box[3]

            # calculate force
            σ = σ₁ + σ₂
            ϵ = ϵ₁ * ϵ₂

            e, ∂e∂x, ∂e∂y, ∂e∂z = lennard_jones_force(vx, vy, vz, σ, ϵ)

            # apply cutoff using a mask for SIMD
            r² = vx * vx + vy * vy + vz * vz
            mask = r² < rcut²

            e *= mask
            ∂e∂x *= mask
            ∂e∂y *= mask
            ∂e∂z *= mask

            # apply switching function
            r = sqrt(r²)
            s, dsdr = shift(r, rcut)
            ∂e∂x = (∂e∂x * s) + (e * dsdr * vx / r)
            ∂e∂y = (∂e∂y * s) + (e * dsdr * vy / r)
            ∂e∂z = (∂e∂z * s) + (e * dsdr * vz / r)
            e *= s

            # accumulate
            e_thread += e
            fx_thread += ∂e∂x
            fy_thread += ∂e∂y
            fz_thread += ∂e∂z
            fx[j] -= ∂e∂x
            fy[j] -= ∂e∂y
            fz[j] -= ∂e∂z
        end

        # accumulate
        fx[i] += fx_thread
        fy[i] += fy_thread
        fz[i] += fz_thread
        e_threads[Threads.threadid()] += e_thread
    end

    e_sum = sum(e_threads)

    # update forces
    forces = force(system)
    @batch for i in eachindex(forces)
        forces[i] += SVector(fx[i], fy[i], fz[i])
    end
    fill!(system.soa.fx, zero(eltype(system.soa.fx)))
    fill!(system.soa.fy, zero(eltype(system.soa.fy)))
    fill!(system.soa.fz, zero(eltype(system.soa.fz)))

    return e_sum
end

Base.@kwdef struct LennardJonesExceptionForce <: AbstractForce
    indices::Vector{Tuple{Int, Int}} = []
    epsilon::Vector{Float64} = []
    sigma::Vector{Float64} = []
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

function lennard_jones_force(vx, vy, vz, σ, ϵ)
    r² = vx * vx + vy * vy + vz * vz

    σ² = σ^2
    σ² = σ² / r²
    σ⁶ = σ² * σ² * σ²

    e = ϵ * σ⁶ * (σ⁶ - 1.0)
    ∂e∂v = -ϵ * σ⁶ * (12 * σ⁶ - 6.0) / r²
    ∂e∂x = ∂e∂v * vx
    ∂e∂y = ∂e∂v * vy
    ∂e∂z = ∂e∂v * vz

    return e, ∂e∂x, ∂e∂y, ∂e∂z
end
