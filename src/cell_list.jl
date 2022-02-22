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
    D = length(box)
    cells = Vector{CartesianIndex{D}}(undef, natoms)
    list = Vector{Int}(undef, natoms)
    cellsize = rsch * ratio
    ncells = map(b -> floor(Int, b / cellsize), Tuple(box))
    nmax = ceil(Int, 1 / ratio)

    if any(ncells .< 2 * nmax + 1)
        error("Box is too small for the cutoff")
    end

    head = OffsetArray(zeros(Int, ncells), CartesianIndex(ntuple(i -> 0, D)...):CartesianIndex((ncells .- 1)...))
    offsets = CartesianIndices((fill(-nmax : nmax, D)...,))[1 : ((2 * nmax + 1)^D + 1) ÷ 2]
    offsets = filter(offset -> norm(Tuple(offset) .- sign.(Tuple(offset))) < 1 / ratio, offsets)
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
    list::Vector{Vector{Int}}
    n::Vector{Int}
end

function NeighborList(natoms::Int, maxnb::Int)
    list = [zeros(Int, maxnb) for _ in 1 : natoms]
    n = zeros(Int, natoms)
    NeighborList(list, n)
end

function update!(nblist::NeighborList, positions, box, rsch, exclusion, cl::LinkedCellList)
    ncells = size(cl.head)
    rsch² = rsch^2

    fill!(nblist.n, zero(eltype(nblist.n)))

    Threads.@threads for ci in CartesianIndices(cl.head)
        i = cl.head[ci]

        @inbounds while i != 0
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

                    # skip to next j-atom if r > rcut
                    r² = v ⋅ v
                    if r² > rsch²
                        j = cl.list[j]
                        continue
                    end

                    nblist.n[i] += 1
                    nblist.list[i][nblist.n[i]] = j

                    # Next j-atom
                    j = cl.list[j]
                end
            end

            # # sort i-atom's neighbors
            # sort!(view(nblist.list[i], 1 : nblist.n[i]))

            # Next i-atom
            i = cl.list[i]
        end
    end # End loop over all cells
    return nblist
end

function update!(nblist::NeighborList, positions, box, rsch, exclusion, cl::NullCellList)
    natoms = length(positions)

    fill!(nblist.n, zero(eltype(nblist.n)))

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

            nblist.n[i] += 1
            nblist.list[i][nblist.n[i]] = j
        end
    end
    return nblist
end
