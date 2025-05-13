struct VectorRadiance{T, U <: AbstractMatrix{T}} <: AbstractMatrix{T}

    S::U

    function VectorRadiance(T::Type, N::Int)

        # We hard-code to 3 Stokes components (I, Q, U)
        return new{T, Matrix{T}}(
            zeros(T, N, 3)
        )

    end

end

@forward VectorRadiance.S (
    Base.length, Base.iterate, Base.size,
    Base.getindex, Base.setindex!, Base.axes,
    Base.copy, Base.copyto!, Base.broadcast,
    Base.broadcast!, Base.dotview, Base.materialize!,
    Base.broadcasted, Base.zero
    )

struct ScalarRadiance{T, U <: AbstractVector{T}} <: AbstractVector{T}

    S::U

    function ScalarRadiance(T::Type, N::Int)
        return new{T, Vector{T}}(zeros(T, N))
    end

end

@forward ScalarRadiance.S (
    Base.length, Base.iterate, Base.size,
    Base.getindex, Base.setindex!, Base.axes,
    Base.copy, Base.copyto!, Base.broadcast,
    Base.broadcast!, Base.dotview, Base.materialize!,
    Base.broadcasted, Base.zero
    )

Radiance = Union{ScalarRadiance, VectorRadiance}


"""
    Overloaded `getproperty`, such that users can access the data inside a
    `ScalarRadiance` object with `.I`, representing the intensity.
"""
function Base.getproperty(obj::ScalarRadiance, sym::Symbol)
    sym === :I && return Base.getfield(obj, :S)
    return Base.getfield(obj, sym)
end

"""
    Overloaded `getproperty`, such that users can access the data inside a
    `ScalarRadiance` object with `.I`, representing the intensity.
"""
function Base.getproperty(obj::VectorRadiance, sym::Symbol)
    sym === :I && return view(obj, :, 1)
    sym === :Q && return view(obj, :, 2)
    sym === :U && return view(obj, :, 3)
    return Base.getfield(obj, sym)
end


function Base.iterate(obj::ScalarRadiance, state::Int=1)
    state == 1 && return Base.getproperty(obj, :I), state + 1
    state == 2 && nothing
end

function Base.iterate(obj::VectorRadiance, state::Int=1)
    state == 1 && return Base.getproperty(obj, :I), state+1
    state == 2 && return Base.getproperty(obj, :Q), state+1
    state == 3 && return Base.getproperty(obj, :U), state+1
    state == 4 &&  nothing
end

struct StokesIterator{T <: Tuple{Vararg{Radiance}}}

    rads::T

    function StokesIterator(objs::Radiance...)
        T = Tuple{Vararg{Radiance}}
        return new{T}(objs)
    end

end

function Base.iterate(SI::StokesIterator, state::Int=1)
    # Intensity is always iterated over
    state == 1 && return (view(x,:,1) for x in SI.rads), state+1
    # Q, U will only be iterated over if all of SI are VectorRadiances
    any(x isa ScalarRadiance for x in SI.rads) && return nothing
    state == 2 && return (view(x,:,2) for x in SI.rads), state+1
    state == 3 && return (view(x,:,3) for x in SI.rads), state+1
    state == 4 && return nothing
end


