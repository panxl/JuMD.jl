const SQRTPI = sqrt(pi)

struct Ewald{D, N}
    alpha::Float64
    kvectors::Vector{NTuple{D, Int}}
    kfactors::Vector{Float64}
    eir::OffsetArray{ComplexF64, 3, Array{ComplexF64, 3}}
    _qeir::NTuple{N, Vector{ComplexF64}}
end

function Ewald(alpha, natoms, kmax, box)
    D = length(box)
    kvectors = Tuple.(CartesianIndices((fill(-kmax : kmax, D)...,))[((2 * kmax + 1)^D + 3) ÷ 2 : end])
    kvectors = filter(x -> (x ⋅ x) <= kmax^2, kvectors)

    kfactors = Vector{Float64}(undef, length(kvectors))

    b = -1 / (4 * alpha^2)
    recip_box = 2 * pi ./ box

    for i in eachindex(kvectors, kfactors)
        k = kvectors[i]
        kr = k .* recip_box
        kr² = kr ⋅ kr
        kfactors[i] = 4 * pi * exp(b * kr²) / kr²
    end

    eir = zeros(ComplexF64, 1:natoms, -kmax:kmax, 1:D)
    _qeir = ntuple(i -> zeros(ComplexF64, natoms), Threads.nthreads())
    Ewald(alpha, kvectors, kfactors, eir, _qeir)
end

function update!(ewald::Ewald, r_box)
    eir = ewald.eir
    kmax = (size(eir)[2] - 1) ÷ 2
    natoms = size(eir)[1]

    # k = 0 and 1
    for i in 1 : natoms
        eir[i, 0, :] .= one(ComplexF64)
        eir[i, 1, :] = cos.(2 * pi * r_box[i]) + sin.(2 * pi * r_box[i])im
    end

    # k from 2 to kmax
    for m in 1 : 3
        for j in 2 : kmax
            for i in 1 : natoms
                eir[i, j, m] = eir[i, j - 1, m] * eir[i, 1, m]
            end
        end
    end

    # negative k values are complex conjugates of positive ones
    for m in 1 : 2
        for j in 1 : kmax
            for i in 1 : natoms
                eir[i, -j, m] = conj(eir[i, j, m])
            end
        end
    end
end

function update!(ewald::Ewald, r_box, box)
    update!(ewald::Ewald, r_box)

    alpha = ewald.alpha
    kvectors = ewald.kvectors
    kfactors = ewald.kfactors

    b = -1 / (4 * alpha^2)
    recip_box = 2 * pi ./ box

    for i in eachindex(kvectors, kfactors)
        k = kvectors[i]
        kr = k .* recip_box
        kr² = kr ⋅ kr
        kfactors[i] = 4 * pi * exp(b * kr²) / kr²
    end
end

function ewald_real_potential(v, α)
    r² = v ⋅ v
    r = sqrt(r²)
    e =  erfc(α * r) / r
    ∂e∂v = -(e + 2 * α * exp(-(α * r)^2) / SQRTPI) / r² * v
    return e, ∂e∂v
end

function ewald_reciporical_potential(v, α)
    r² = v ⋅ v
    r = sqrt(r²)
    e = erf(α * r) / r
    ∂e∂v = -(e - 2 * α * exp(-(α * r)^2) / SQRTPI) / r² * v
    return e, ∂e∂v
end
