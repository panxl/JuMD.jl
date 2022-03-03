abstract type AbstractCellList end

struct NullCellList <: AbstractCellList end

function update!(::NullCellList, positions, box) end

struct LinkedCellList{D} <: AbstractCellList
    cells::Vector{CartesianIndex{D}}
    list::Vector{Int}
    head::OffsetArray{Int, D, Array{Int, D}}
    offsets::Vector{CartesianIndex{D}}
end

function LinkedCellList(natoms, rsch, box; ratio=1.0)
    ndim = length(box)
    cells = Vector{CartesianIndex{ndim}}(undef, natoms)
    list = Vector{Int}(undef, natoms)
    ncells = map(b -> floor(Int, b / (rsch * ratio)), Tuple(box))
    ncells = floor.(Int, ncells ./ 2) .* 2
    cellsize = box ./ ncells
    nmax = ceil.(Int, rsch ./ cellsize)

    if any(ncells .< 2 .* nmax .+ 1)
        error("Box is too small for the cutoff")
    end

    head = OffsetArray(zeros(Int, ncells), CartesianIndex(ntuple(i -> 0, ndim)...):CartesianIndex((ncells .- 1)...))
    offsets = CartesianIndices(Tuple(fill(-n : n)[] for n in nmax))[1 : (prod(2 .* nmax .+ 1) + 1) ÷ 2]
    offsets = filter(offset -> norm((Tuple(offset) .- sign.(Tuple(offset))) .* cellsize) < rsch, offsets)
    LinkedCellList(cells, list, head, offsets)
end

function update!(cl::LinkedCellList, positions, box)
    inv_box = 1 ./ box
    ncells = size(cl.head)
    fill!(cl.head, zero(eltype(cl.head)))
    @inbounds for i in eachindex(positions)
        r_box = positions[i] .* inv_box
        ci = find_cell(r_box, ncells)
        cl.cells[i] = ci
        cl.list[i] = cl.head[ci]
        cl.head[ci] = i
    end
    return cl
end

# find cell index given fractional coordinates
function find_cell(r_box, ncells)
    r_box = wrap_positions(r_box)
    cell = CartesianIndex(Tuple(map(x -> trunc(Int, x), r_box .* ncells)))
    return cell
end

# put fractional coordinates in the range of [0, 1]
wrap_positions(r_box) = r_box .- floor.(r_box)

# put fractional coordinates in the range of [-0.5, 0.5]
minimum_image(r_box) = r_box .- round.(r_box)

struct NeighborList
    n_neighbors::Vector{Int}
    lists::Vector{Vector{Int}}
    rskin::Float64
    x0::Vector{Float64}
    y0::Vector{Float64}
    z0::Vector{Float64}
end

function NeighborList(natoms::Int, maxnb::Int; rskin=0.0)
    n_neighbors = zeros(Int, natoms)
    lists = [zeros(Int, maxnb) for _ in 1 : natoms]
    x0 = zeros(Float64, natoms)
    y0 = zeros(Float64, natoms)
    z0 = zeros(Float64, natoms)
    NeighborList(n_neighbors, lists, rskin, x0, y0, z0)
end

function update!(nbl::NeighborList, positions, box, rcut, exclusion, x, y, z, cl::LinkedCellList)
    # check if update is needed
    threshold = (0.5 * nbl.rskin)^2
    @unpack x0, y0, z0 = nbl
    n = 0
    @tturbo for i in eachindex(x, y, z, x0, y0, z0)
        d² = (x[i] - x0[i])^2 + (y[i] - y0[i])^2 + (z[i] - z0[i])^2
        n += d² > threshold
    end

    # skip updating neighbor list because it is not needed
    if n == 0
        return nbl
    end

    # update saved positions
    @tturbo for i in eachindex(x, y, z, x0, y0, z0)
        x0[i] = x[i]
        y0[i] = y[i]
        z0[i] = z[i]
    end

    # now we are updating cell list and neighbor list
    update!(cl, positions, box)

    @unpack n_neighbors, lists, rskin = nbl

    ncells = size(cl.head)
    rsch = rcut + rskin
    rsch² = rsch^2

    # empty the current neighbor list
    fill!(n_neighbors, zero(eltype(n_neighbors)))

    @batch for ci in CartesianIndices(cl.head)
        i = cl.head[ci]

        while i != 0
            x1 = positions[i]

            for offset in cl.offsets
                if offset == last(cl.offsets)
                    j = cl.list[i]
                else
                    cj = Tuple(ci + offset)
                    cj = mod.(cj, ncells)
                    j = cl.head[CartesianIndex(cj)]
                end

                while j != 0
                    # skip if j-atom is in i-atom's exclusion list
                    if !isnothing(exclusion) && j ∈ exclusion[i]
                        j = cl.list[j]
                        continue
                    end

                    x2 = positions[j]
                    v = x2 - x1

                    # apply minimum image convention
                    if any(box .!= 0.0)
                        v = minimum_image(v ./ box) .* box
                    end

                    # skip to next j-atom if r > rsch
                    r² = v ⋅ v
                    if r² > rsch²
                        j = cl.list[j]
                        continue
                    end

                    # add j-atom to i-atom's neighbor list
                    n_neighbors[i] += 1
                    lists[i][n_neighbors[i]] = j

                    # Next j-atom
                    j = cl.list[j]
                end
            end

            # Next i-atom
            i = cl.list[i]
        end
    end # End loop over all cells
    return nbl
end

function update!(nbl::NeighborList, positions, box, rsch, exclusion, x, y, z, cl::NullCellList)
    @unpack n_neighbors, lists = nbl
    natoms = length(positions)

    # empty the current neighbor list
    fill!(n_neighbors, zero(eltype(n_neighbors)))

    for i in 1 : (natoms - 1)
        x1 = positions[i]

        for j in (i + 1 : natoms)
            # skip if j-atom is in i-atom's exclusion list
            if j ∈ exclusion[i]
                continue
            end

            x2 = positions[j]
            v = x2 - x1

            # apply minimum image convention
            if any(box .!= 0.0)
                v = minimum_image(v ./ box) .* box
            end

            # skip to next j-atom if r > rsch
            if !isnothing(rsch)
                r = norm(v)
                if r > rsch
                    continue
                end
            end

            n_neighbors[i] += 1
            lists[i][n_neighbors[i]] = j
        end
    end
    return nbl
end
