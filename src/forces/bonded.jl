Base.@kwdef struct HarmonicBondForce <: AbstractForce
    indices::Vector{Tuple{Int, Int}} = []
    length::Vector{Float64} = []
    k::Vector{Float64} = []
end

function force!(system, f::HarmonicBondForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.length, f.k)
        i1, i2 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        r₀= f.length[i]
        k = f.k[i]

        (f1, f2), e = harmonic_bond_force(x1, x2, r₀, k)

        forces[i1] += f1
        forces[i2] += f2
        e_sum += e
    end
    return e_sum
end

function harmonic_bond_force(x1, x2, r₀, k)
    v = x2 - x1
    r = norm(v)
    Δr = r - r₀
    f1 = k * Δr / r * v
    f2 = -f1
    e = 0.5 * k * Δr^2
    return (f1, f2), e
end

Base.@kwdef struct HarmonicAngleForce <: AbstractForce
    indices::Vector{Tuple{Int, Int, Int}} = []
    angle::Vector{Float64} = []
    k::Vector{Float64} = []
end

function force!(system, f::HarmonicAngleForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.angle, f.k)
        i1, i2, i3 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        x3 = positions[i3]
        θ₀= f.angle[i]
        k = f.k[i]

        (f1, f2, f3), e = harmonic_angle_force(x1, x2, x3, θ₀, k)

        forces[i1] += f1
        forces[i2] += f2
        forces[i3] += f3
        e_sum += e
    end
    return e_sum
end

function harmonic_angle_force(x1, x2, x3, θ₀, k)
    v1 = x2 - x1
    v2 = x2 - x3
    r1 = norm(v1)
    r2 = norm(v2)
    cosine = clamp((v1 ⋅ v2) / (r1 * r2), -1, 1)
    θ = acos(cosine)
    Δθ = θ - θ₀
    dEdθ = k * Δθ

    cp = v1 × v2
    rp = norm(cp)
    f1 = dEdθ / (r1 * r1 * rp) * (v1 × cp)
    f3 = dEdθ / (r2 * r2 * rp) * (cp × v2)
    f2 = -(f1 + f3)
    e = 0.5 * k * Δθ^2
    return (f1, f2, f3), e
end

Base.@kwdef struct PeriodicTorsionForce <: AbstractForce
    indices::Vector{Tuple{Int, Int, Int, Int}} = []
    periodicity::Vector{Int} = []
    phase::Vector{Float64} = []
    k::Vector{Float64} = []
end

function force!(system, f::PeriodicTorsionForce)
    positions = position(system)
    forces = force(system)
    e_sum = 0.0

    @inbounds for i in eachindex(f.indices, f.periodicity, f.phase, f.k)
        i1, i2, i3, i4 = f.indices[i]

        x1 = positions[i1]
        x2 = positions[i2]
        x3 = positions[i3]
        x4 = positions[i4]
        n = f.periodicity[i]
        ϕ₀ = f.phase[i]
        k = f.k[i]

        (f1, f2, f3, f4), e = periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)

        forces[i1] += f1
        forces[i2] += f2
        forces[i3] += f3
        forces[i4] += f4
        e_sum += e
    end
    return e_sum
end

function periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)
    v1 = x2 - x1
    v2 = x2 - x3
    v3 = x4 - x3

    r₂² = v2 ⋅ v2
    r₂ = sqrt(r₂²)
    cp1 = v1 × v2
    cp2 = v2 × v3
    x = cp1 ⋅ cp2
    y = (cp1 × cp2) ⋅ v2 / r₂
    ϕ = atan(y, x)

    Δϕ = n * ϕ - ϕ₀
    dEdϕ = -k * n * sin(Δϕ)

    f1 =  dEdϕ * r₂ / (cp1 ⋅ cp1) * cp1
    f4 = -dEdϕ * r₂ / (cp2 ⋅ cp2) * cp2

    s = (v1 ⋅ v2) / r₂² * f1 - (v3 ⋅ v2) / r₂² * f4
    f2 =  s - f1
    f3 = -s - f4

    e = k * (1.0 + cos(Δϕ))
    return (f1, f2, f3, f4), e
end
