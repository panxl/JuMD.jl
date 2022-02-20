abstract type AbstractThermostat end

struct NullThermostat <: AbstractThermostat end

function (::NullThermostat)(system, dt) end

struct LangevinThermostat <: AbstractThermostat
    gamma::Float64
    T::Float64
end

function (thermostat::LangevinThermostat)(system, dt)
    kT = KB * thermostat.T
    v = velocity(system)
    inv_M = inverse_mass(system)
    γ = thermostat.gamma

    c1 = exp(-γ * dt)
    c2 = sqrt(1 - c1^2)
    @inbounds for i in eachindex(v, inv_M)
        v[i] = c1 * v[i] + c2 * sqrt(kT * inv_M[i]) * randn(eltype(v))
    end
    return nothing
end

struct GJFThermostat <: AbstractThermostat
    gamma::Float64
    T::Float64
end

function (thermostat::GJFThermostat)(system, dt)
    kT = KB * thermostat.T
    v = velocity(system)
    inv_M = inverse_mass(system)
    γ = thermostat.gamma

    c1 = (1 - 0.5 * γ * dt) / (1 + 0.5 * γ * dt)
    c2 = sqrt(1 - c1^2)
    @inbounds for i in eachindex(v, inv_M)
        v[i] = c1 * v[i] + c2 * sqrt(kT * inv_M[i]) * randn(eltype(v))
    end
    return nothing
end

function random_velocities!(v, inv_M, kT)
    @inbounds for i in eachindex(v, inv_M)
        v[i] = sqrt(kT * inv_M[i]) * randn(eltype(v))
    end
    return nothing
end