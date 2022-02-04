abstract type AbstractReporter end

struct StateReporter <: AbstractReporter
    interval::Int
    path::String
end

function (reporter::StateReporter)(io::IO, system::AbstractSystem, current_step::Integer)
    if current_step == 1
        println(io, "# Step    T    Etot    Epot    Ekin")
    end
    if current_step % reporter.interval == 0
        masses = mass(system)
        velocities = velocity(system)
        velocities_half = velocity_half(system)
        ekin = 0.0
        for (m, v) in zip(masses, velocities)
            ekin += 0.5 * m * (v ⋅ v)
        end
        ekin_half = 0.0
        for (m, v) in zip(masses, velocities_half)
            ekin_half += 0.5 * m * (v ⋅ v)
        end
        T = 2 * ekin_half / (3 * length(masses) * KB)
        epot = sum(system.force_groups.energies)
        etot = epot + ekin
        println(io, current_step, " ", T, " ",  etot, " ", epot, " ", ekin)
        return nothing
    end
end

struct CheckPointReporter <: AbstractReporter
    interval::Int
end

struct XYZTrajectoryReporter <: AbstractReporter
    interval::Int
    path::String
end

function (reporter::XYZTrajectoryReporter)(io::IO, system::AbstractSystem, current_step)
    if current_step % reporter.interval == 0
        positions = position(system)
        atomic_numbers = atomic_number(system)
        println(io, length(positions), "\n")
        for (n, x) in zip(atomic_numbers, positions)
            println(io, n, " ", x[1], " ", x[2], " ", x[3])
        end
    end
end

struct CVReporter <: AbstractReporter
    interval::Int
end