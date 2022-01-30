import LinearAlgebra: norm, cross, dot

abstract type AbstractSystem{D} end

abstract type AbstractForce end

struct HarmonicBondForce <: AbstractForce
    indices::Tuple{Int, Int}
    r0::Float64
    k::Float64
end

function force!(system::AbstractSystem, f::HarmonicBondForce)
    i1, i2 = f.indices
    x1 = position(system)[i1]
    x2 = position(system)[i2]
    r0= f.r0
    k = f.k

    f1, f2, e = harmonic_bond_force(x1, x2, r0, k)
    force(system)[i1] += f1
    force(system)[i2] += f2

    return e
end

function harmonic_bond_force(x1, x2, r0, k)
    v21 = x2 - x1
    r = norm(v21)
    Δr = r - r0
    f1 = k * Δr / r * v21
    f2 = -f1
    e = 0.5 * k * Δr^2
    return f1, f2, e
end

struct HarmonicAngleForce <: AbstractForce
    indices::Tuple{Int, Int, Int}
    a0::Float64
    k::Float64
end

function force!(system::AbstractSystem, f::HarmonicAngleForce)
    i1, i2, i3 = f.indices
    x1 = position(system)[i1]
    x2 = position(system)[i2]
    x3 = position(system)[i3]
    a0= f.a0
    k = f.k

    f1, f2, f3, e = harmonic_angle_force(x1, x2, x3, a0, k)
    force(system)[i1] += f1
    force(system)[i2] += f2
    force(system)[i3] += f3

    return e
end

function harmonic_angle_force(x1, x2, x3, a0, k)
    v21 = x2 - x1
    v23 = x2 - x3
    r21 = norm(v21)
    r23 = norm(v23)
    cosine = dot(v21, v23) / (r21 * r23)
    a = acos(cosine)
    Δa = a - a0

    cp = cross(v21, v23)
    rcp = norm(cp)
    f1 = k * Δa / (r21 * r21 * rcp) * cross(v21, cp)
    f3 = k * Δa / (r23 * r23 * rcp) * cross(cp, v23)
    f2 = -(f1 + f3)
    e = 0.5 * k * Δa^2
    return f1, f2, f3, e
end