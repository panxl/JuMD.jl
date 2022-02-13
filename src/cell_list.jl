abstract type AbstractCellList end

struct NullCellList <: AbstractCellList end

struct LinkedCellList{D} <: AbstractCellList
    cells::Vector{CartesianIndex{D}}
    list::Vector{Int}
    head::OffsetArray{Int, D, Array{Int, D}}
    offsets::Vector{CartesianIndex{D}}
    cutoff::Float64
end

function LinkedCellList(natoms, cutoff, box)
    D = length(box)
    cells = Vector{CartesianIndex{D}}(undef, natoms)
    list = Vector{Int}(undef, natoms)
    ncells = map(b -> floor(Int, b / cutoff), box)
    head = OffsetArray(zeros(Int, ncells), CartesianIndex(ntuple(i->0, D)...):CartesianIndex((ncells .- 1)...))
    offsets = CartesianIndices((fill(-1:1, D)...,))[1:(3^D + 1)÷2]
    LinkedCellList(cells, list, head, offsets, cutoff)
end

function update!(cl::LinkedCellList, positions, box)
    inv_box = 1 ./ box
    ncells = size(cl.head)
    fill!(cl.head, zero(eltype(cl.head)))
    @inbounds for i in eachindex(positions)
        ci = find_cell((positions[i] .* inv_box), ncells)
        cl.cells[i] = ci
        cl.list[i] = cl.head[ci]
        cl.head[ci] = i
    end
    return cl
end

function update!(::NullCellList, positions, box) end

function find_cell(r_box, ncells)
    r_box = wrap_positions(r_box)
    cell = CartesianIndex(Tuple(map(x -> trunc(Int, x), r_box .* ncells)))
    return cell
end

wrap_positions(r_box) = r_box .- floor.(r_box)

minimum_image(r_box) = r_box .- round.(r_box)
