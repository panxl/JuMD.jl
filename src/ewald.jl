abstract type AbstractRecip end

struct NullRecip <: AbstractRecip end

struct EwaldRecip{D, N} <: AbstractRecip
    alpha::Float64
    kvectors::Vector{NTuple{D, Int}}
    eir::OffsetArray{ComplexF64, 3, Array{ComplexF64, 3}}
    _qeir::NTuple{N, Vector{ComplexF64}}
end

function EwaldRecip(alpha, natoms, kmax, box)
    D = length(box)
    kvectors = Tuple.(CartesianIndices((fill(-kmax : kmax, D)...,))[((2 * kmax + 1)^D + 3) ÷ 2 : end])
    kvectors = filter(x -> (x ⋅ x) <= kmax^2, kvectors)
    eir = zeros(ComplexF64, 1:natoms, -kmax:kmax, 1:D)
    _qeir = ntuple(i -> zeros(ComplexF64, natoms), Threads.nthreads())
    EwaldRecip(alpha, kvectors, eir, _qeir)
end

function EwaldRecip(box, cutoff, natoms; tol=1e-5, alpha=nothing, kmax=nothing)
    if isnothing(alpha)
        alpha = get_ewald_alpha(cutoff, tol)
    end
    if isnothing(kmax)
        # hard code for now
        kmax = 20
    end
    EwaldRecip(alpha, natoms, kmax, box)
end

function update!(recip::EwaldRecip, positions, box)
    eir = recip.eir
    kmax = (size(eir)[2] - 1) ÷ 2
    natoms = size(eir)[1]

    @inbounds for m in 1 : 3
        # k = 0 and 1
        for i in 1 : natoms
            r_box = positions[i][m] / box[m]
            r_box = 2 * pi * (r_box - round(r_box))

            eir[i, 0, m] = one(ComplexF64)
            eir[i, 1, m] = cos(r_box) + sin(r_box)im
        end

        # k from 2 to kmax
        for j in 2 : kmax
            for i in 1 : natoms
                eir[i, j, m] = eir[i, j - 1, m] * eir[i, 1, m]
            end
        end
    end

    # negative k values are complex conjugates of positive ones
    @inbounds for m in 1 : 2
        for j in 1 : kmax
            for i in 1 : natoms
                eir[i, -j, m] = conj(eir[i, j, m])
            end
        end
    end
end

function (recip::EwaldRecip)(box, charges, positions, force_array)
    update!(recip, positions, box)

    alpha = recip.alpha
    kvectors = recip.kvectors
    eir = recip.eir

    b = -1 / (4 * alpha^2)
    recip_box = 2 * pi ./ box
    inv_v = 1 / prod(box)

    e_threads = zeros(Threads.nthreads())

    Threads.@threads for i in eachindex(recip.kvectors)
        forces = force_array[Threads.threadid()]
        qeir = recip._qeir[Threads.threadid()]

        k = kvectors[i]
        kx, ky, kz = k

        kr = k .* recip_box
        kr² = kr ⋅ kr
        kfac = 4 * pi * exp(b * kr²) / kr²

        term = zero(ComplexF64)

        @inbounds @simd for j in 1 : length(charges)
            qeir[j] = charges[j] * eir[j, kx, 1] * eir[j, ky, 2] * eir[j, kz, 3]
            term += qeir[j]
        end

        e_threads[Threads.threadid()] += kfac * real(conj(term) * term)

        @inbounds @simd for j in 1 : length(charges)
            f_recip = 2 * kfac * imag(conj(term) * qeir[j]) .* kr
            forces[j] += f_recip .* (KE * inv_v)
        end
    end

    e_sum = sum(e_threads) * (KE * inv_v)

    # substract self part of k-space sum
    e_sum -= KE * alpha * (charges ⋅ charges) / SQRTPI

    return e_sum
end

function get_ewald_alpha(cutoff, tol)
    alpha = 1
    while erfc(alpha * cutoff) / cutoff >= tol
        alpha *= 2
    end
    alpha_lo = 0
    alpha_hi = alpha
    for _ in 1:100
        alpha = .5 * (alpha_lo + alpha_hi)
        if erfc(alpha * cutoff) / cutoff >= tol
            alpha_lo = alpha
        else
            alpha_hi = alpha
        end
    end
    return alpha
end
