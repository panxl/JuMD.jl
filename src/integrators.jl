abstract type AbstractIntegrator end

struct VerlocityVerletIntegrator{T,C} <: AbstractIntegrator
    dt::Float64
    thermostat::T
    constraint::C
end

function (integrator::VerlocityVerletIntegrator)(system::AbstractSystem)
    x = position(system)
    v = velocity(system)
    F = force(system)
    inv_M = inverse_mass(system)
    halfdt = integrator.dt / 2

    @. v = v + F * inv_M * halfdt
    integrator.constraint(system)
    @. x = x + v * halfdt
    integrator.thermostat(system, integrator.dt)
    @. x = x + v * halfdt
    integrator.constraint(system)
    force!(system)
    integrator.constraint(system)
    @. v = v + F * inv_M * halfdt
    integrator.constraint(system)
    return nothing
end