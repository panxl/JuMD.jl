include("helpme/helpme.jl")
using .helpme

struct PMERecip <: AbstractRecip
    pme::Ptr{PMEInstance}
    alpha::Float64
    spline_order::Int
    grid::Tuple{Int, Int, Int}
    scaling::Float64
    nthreads::Int
end

function PMERecip(alpha, spline_order, grid, scaling, nthreads)
    pme = helpme_createD()
    helpme_setupD(pme, 1, alpha, spline_order, grid..., scaling, nthreads)
    PMERecip(pme, alpha, spline_order, grid, scaling, nthreads)
end

function update!(recip::PMERecip, box)
    helpme_set_lattice_vectorsD(recip.pme, box..., 90, 90, 90, XAligned)
end

function (recip::PMERecip)(box, charges, positions, force_array)
    update!(recip, box)

    natoms = length(charges)
    e = helpme_compute_EF_recD(recip.pme, natoms, 0, charges, reinterpret(Float64, positions), reinterpret(Float64, force_array[1]))

    # substract self part of k-space sum
    e -= recip.scaling * recip.alpha * (charges ⋅ charges) / SQRTPI

    return e
end

function potential!(recip::PMERecip, charges, positions, sites, potentials)
    natoms = length(charges)
    nsites = length(sites)
    helpme_compute_P_recD(recip.pme, natoms, 0, charges, reinterpret(Float64, positions), nsites, reinterpret(Float64, sites), 1, reinterpret(Float64, potentials))
end
