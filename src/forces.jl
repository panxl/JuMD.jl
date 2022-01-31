abstract type AbstractSystem{D} end

abstract type AbstractForce end

struct HarmonicBondForce <: AbstractForce
    indices::Tuple{Int, Int}
    length::Float64
    k::Float64
end

function force!(system::AbstractSystem, f::HarmonicBondForce)
    i1, i2 = f.indices
    x1 = position(system)[i1]
    x2 = position(system)[i2]
    r₀= f.length
    k = f.k

    (f1, f2), e = harmonic_bond_force(x1, x2, r₀, k)
    force(system)[i1] += f1
    force(system)[i2] += f2

    return e
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

struct HarmonicAngleForce <: AbstractForce
    indices::Tuple{Int, Int, Int}
    angle::Float64
    k::Float64
end

function force!(system::AbstractSystem, f::HarmonicAngleForce)
    i1, i2, i3 = f.indices
    x1 = position(system)[i1]
    x2 = position(system)[i2]
    x3 = position(system)[i3]
    θ₀= f.angle
    k = f.k

    (f1, f2, f3), e = harmonic_angle_force(x1, x2, x3, θ₀, k)
    force(system)[i1] += f1
    force(system)[i2] += f2
    force(system)[i3] += f3

    return e
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

struct PeriodicTorsionForce <: AbstractForce
    indices::Tuple{Int, Int, Int, Int}
    periodicity::Int
    phase::Float64
    k::Float64
end

function force!(system::AbstractSystem, f::PeriodicTorsionForce)
    i1, i2, i3, i4 = f.indices
    x1 = position(system)[i1]
    x2 = position(system)[i2]
    x3 = position(system)[i3]
    x4 = position(system)[i4]
    n = f.periodicity
    ϕ₀ = f.phase
    k = f.k

    (f1, f2, f3, f4), e = periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)
    force(system)[i1] += f1
    force(system)[i2] += f2
    force(system)[i3] += f3
    force(system)[i4] += f4

    return e
end

function periodic_torsion_force(x1, x2, x3, x4, n, ϕ₀, k)
    v1 = x2 - x1
    v2 = x2 - x3
    v2 /= norm(v2)
    v3 = x4 - x3
    cp1 = v1 × v2
    cp2 = v2 × v3
    m1 = cp1 × v2
    x = cp1 ⋅ cp2
    y = m1 ⋅ cp2
    ϕ = atan(y, x)

    Δϕ = n * ϕ - ϕ₀
    dEdϕ = -k * n * sin(Δϕ)

    f1 = -dEdϕ / (cp1 ⋅ cp1) * cp1
    f4 = dEdϕ / (cp2 ⋅ cp2) * cp2

    s = (v1 ⋅ v2) * f1 - (v3 ⋅ v2) * f4
    f2 = s - f1
    f3 = -s - f4

    e = k * (1.0 + cos(Δϕ))
    return (f1, f2, f3, f4), e
end