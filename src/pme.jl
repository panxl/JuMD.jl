include("helpme/helpme.jl")
using .helpme

struct PMERecip <: AbstractRecip
    pme::Ptr{PMEInstance}
    alpha::Float64
    nfft::Tuple{Int, Int, Int}
    spline_order::Int
    scaling::Float64
    nthreads::Int
end

function PMERecip(alpha, nfft, spline_order, scaling, nthreads)
    pme = helpme_createD()
    helpme_setupD(pme, 1, alpha, spline_order, nfft..., scaling, nthreads)
    PMERecip(pme, alpha, nfft, spline_order, scaling, nthreads)
end

function PMERecip(box, cutoff; tol=1e-5, alpha=nothing, grid_spacing=0.1, spline_order=6, scaling=1.0, nthreads=1)
    if isnothing(alpha)
        alpha = get_ewald_alpha(cutoff, tol)
    end
    nfft = Tuple(find_fft_dimension.(ceil.(Int, box ./ grid_spacing)))
    PMERecip(alpha, nfft, spline_order, scaling, nthreads)
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

function find_fft_dimension(minimum)
    if minimum < 1
        return 1
    end
    while true
        unfactored = minimum
        for factor in 2 : 5
            while unfactored > 1 && unfactored % factor == 0
                unfactored /= factor
                if unfactored == 1
                    return minimum
                end
            end
        end
    minimum += 1
    end
end
