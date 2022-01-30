abstract type AbstractReporter end

struct StateReporter{T<:Integer} <: AbstractReporter
    interval::T
end

struct CheckPointReporter{T<:Integer} <: AbstractReporter
    interval::T
end

struct NetCDFReporter{T<:Integer} <: AbstractReporter
    interval::T
end

struct CVReporter{T<:Integer} <: AbstractReporter
    interval::T
end