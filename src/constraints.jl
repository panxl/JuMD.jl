abstract type AbstractConstraint end

struct NullConstraint <: AbstractConstraint end

function (::NullConstraint)(system::AbstractSystem) end

struct SHAKE <: AbstractConstraint
    tol::Float64
end

struct RATTLE <: AbstractConstraint
    tol::Float64
end