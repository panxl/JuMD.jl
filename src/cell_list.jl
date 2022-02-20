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
