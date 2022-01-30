abstract type AbstractSimulation end

struct Simulation{S<:AbstractSystem,I<:AbstractIntegrator,R<:AbstractReporter} <: AbstractSimulation
    system::S
    integrator::I
    reporters::Vector{R}
end

function Base.show(io::IO, ::MIME"text/plain", sim::Simulation)
    # println(io, "System: $(sim.system)")
    println(io, "Integrator: $(sim.integrator)")
    println(io, "Reporters: $(sim.reporters)")
end

function Base.run(sim::AbstractSimulation, steps::Integer)
    for _ in 1:steps
        sim.integrator(sim.system)
    end
end