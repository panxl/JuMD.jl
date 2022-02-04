abstract type AbstractSimulation end

struct Simulation{S<:AbstractSystem,I<:AbstractIntegrator,R} <: AbstractSimulation
    system::S
    integrator::I
    reporters::R
end

function Base.show(io::IO, ::MIME"text/plain", sim::Simulation)
    # println(io, "System: $(sim.system)")
    println(io, "Integrator: $(sim.integrator)")
    println(io, "Reporters: $(sim.reporters)")
end

function Base.run(sim::AbstractSimulation, steps::Integer)
    files = map(reporter::AbstractReporter -> open(reporter.path, "w"), sim.reporters)
    try
        for i in 1:steps
            sim.integrator(sim.system)
            foreach((reporter::AbstractReporter, io::IO) -> reporter(io, sim.system, i), sim.reporters, files)
        end
    finally
        close.(files)
    end
end