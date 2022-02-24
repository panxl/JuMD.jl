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
        e = force!(system, forces, neighbor_list(system))
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
    e_threads = zeros(Threads.nthreads())

    @batch for i in 1 : length(positions)
        forces = force(system, Threads.threadid())
        e_thread = 0.0

        # skip to next i-atom if i-atom's epsilon is zero
        ϵ₁ = f.epsilon[i]
        if ϵ₁ == 0.0
            continue
        end

        # load other i-atom's information
        x1 = positions[i]
        σ₁ = f.sigma[i]

        ilist = nbl.list[i]
        for j_idx in 1 : nbl.n[i]
            j = ilist[j_idx]
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
            r = norm(v)
            s, dsdr = shift(r, rcut)
            ∂e∂v = ∂e∂v .* s + e * dsdr / r * v
            e *= s

            forces[i] += ∂e∂v
            forces[j] -= ∂e∂v
            e_thread += e
        end

        e_threads[Threads.threadid()] += e_thread
    end

    e_sum = sum(e_threads)
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
