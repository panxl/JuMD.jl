abstract type AbstractThermostat end

struct NullThermostat <: AbstractThermostat end

function (::NullThermostat)(system::AbstractSystem, dt) end

struct LangevinThermostat <: AbstractThermostat
    gamma::Float64
    T::Float64
end

function (thermostat::LangevinThermostat)(system::AbstractSystem, dt)
    kT = KB * thermostat.T
    v = velocity(system)
    inv_M = inverse_mass(system)
    langevin_thermostat(v, inv_M, dt, thermostat.gamma, kT)
end

function langevin_thermostat(v, inv_M, dt, gamma, kT)
    c1 = exp(-gamma * dt)
    c2 = sqrt(1 - c1^2)
    @inbounds for i in eachindex(v)
        v[i] = c1 * v[i] + c2 * sqrt(kT * inv_M[i]) * randn(eltype(v))
    end
    return nothing
end

function random_velocities!(v, inv_M, kT)
    @inbounds for i in eachindex(v)
        v[i] = sqrt(kT * inv_M[i]) .* randn(eltype(v))
    end
    return nothing
end