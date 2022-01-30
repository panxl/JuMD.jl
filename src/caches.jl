abstract type AbstractCache end

struct Cache{T} <: AbstractCache 
    data::Dict{Symbol, Vector{T}}
end

Cache(T::Type) = Cache(Dict{Symbol, Vector{T}}())

Base.hasproperty(c::Cache, name::Symbol) = haskey(getfield(c, :data), name)
Base.getproperty(c::Cache, name::Symbol) = getfield(c, :data)[name]
Base.setproperty!(c::Cache, name::Symbol, value) = getfield(c, :data)[name] = value
Base.propertynames(c::Cache) = keys(getfield(c, :data))