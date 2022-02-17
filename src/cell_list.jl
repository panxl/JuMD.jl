abstract type AbstractCellList end

struct NullCellList <: AbstractCellList end

struct LinkedCellList{D} <: AbstractCellList
    cells::Vector{CartesianIndex{D}}
    list::Vector{Int}
    head::OffsetArray{Int, D, Array{Int, D}}
    offsets::Vector{CartesianIndex{D}}
end

function LinkedCellList(natoms, cellsize, box)
    D = length(box)
    cells = Vector{CartesianIndex{D}}(undef, natoms)
    list = Vector{Int}(undef, natoms)
    ncells = map(b -> floor(Int, b / cellsize), Tuple(box))

    if any(ncells .< 3)
        error("Box is too small for the cutoff")
    end

    head = OffsetArray(zeros(Int, ncells), CartesianIndex(ntuple(i -> 0, D)...):CartesianIndex((ncells .- 1)...))
    offsets = CartesianIndices((fill(-1 : 1, D)...,))[1 : (3^D + 1) ÷ 2]
    LinkedCellList(cells, list, head, offsets)
end

function update!(cl::LinkedCellList, r_box)
    ncells = size(cl.head)
    fill!(cl.head, zero(eltype(cl.head)))
    @inbounds for i in eachindex(r_box)
        ci = find_cell(r_box[i], ncells)
        cl.cells[i] = ci
        cl.list[i] = cl.head[ci]
        cl.head[ci] = i
    end
    return cl
end

function update!(::NullCellList, r_box) end

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
