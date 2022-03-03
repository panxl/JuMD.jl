abstract type AbstractIntegrator end

Base.Base.@kwdef struct VerlocityVerletIntegrator{T, C} <: AbstractIntegrator
    dt::Float64
    thermostat::T = NullThermostat()
    constraint::C = NullConstraint()
end

function (integrator::VerlocityVerletIntegrator)(system)
    x = position(system)
    v = velocity(system)
    x_last = position_last(system)
    v_half = velocity_half(system)
    F = force(system)
    inv_M = inverse_mass(system)
    halfdt = integrator.dt / 2

    @. v = v + F * inv_M * halfdt
    integrator.constraint(system)

    x_last .= x
    @. x = x + v * halfdt
    integrator.thermostat(system, integrator.dt)
    @. x = x + v * halfdt
    integrator.constraint(system)
    force!(system)
    integrator.constraint(system)

    # save half step velocities
    if integrator.thermostat isa GJFThermostat
        dt = integrator.dt
        γ = integrator.thermostat.gamma
        @. v_half = sqrt(1 + 0.5 * γ * dt) / dt * (x - x_last)
    else
        v_half .= v
    end

    @. v = v + F * inv_M * halfdt
    integrator.constraint(system)

    update!(system.soa, x)

    return nothing
end
