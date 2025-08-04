"""
$(SIGNATURES)

Pretty printing for State Vector types
"""
function show(io::IO, ::MIME"text/plain", SV::RetrievalStateVector)

    pretty_table(
        io,
        permutedims(hcat([
            ["#$(i)", get_name(sve), get_current_value(sve), get_unit(sve)]
              for (i, sve) in enumerate(SV.state_vector_elements)]...
                 )),
        noheader=true, alignment=[:l, :r, :r, :l], title="State Vector (current)",
        hlines=[0,], vlines=[1,], display_size=(150, 100)
    )

    # "old" way
    #for (i, sve) in enumerate(SV.state_vector_elements)
    #    println("#$(i), $(get_name(sve)): $(sve.iterations[end])")
    #end

end

"""
$(TYPEDSIGNATURES)

Basic printing for State Vector types
"""
function show(io::IO, SV::RetrievalStateVector)

    for (i, sve) in enumerate(SV.state_vector_elements)
        println("#$(i), $(get_name(sve)): $(sve.iterations[end])")
    end
end

"""
Gets the names of all state vector elements in the state vector

$(TYPEDSIGNATURES)
"""
get_name(SV::RetrievalStateVector) = map(get_name, SV.state_vector_elements)


"""
$(TYPEDSIGNATURES)

Return the length of the state vector
"""
function length(SV::RetrievalStateVector)

    return length(SV.state_vector_elements)

end



function get_posterior_ucert(q::OEQuantities, SVE::AbstractStateVectorElement)

    for (idx, qSVE) in enumerate(q.state_vector.state_vector_elements)
        if qSVE === SVE
            return q.SV_ucert[idx]
        end
    end

end


"""
$(TYPEDSIGNATURES)

Returns all state vector elements from state vector `SV` that satisfy the type `t`
"""
function return_SVE_of_type(
    SV::RetrievalStateVector,
    t::Type{T}
    ) where T<:AbstractStateVectorElement

    # Allocate result container
    result = AbstractStateVectorElement[]

    # Loop through SVEs that match the type
    for (idx, sve) in StateVectorIterator(SV, t)
        push!(result, sve)
    end

    return result

end


"""
$(TYPEDSIGNATURES)

Returns true if any of the state vector elements
is of the type requested. Otherwise returns false.
"""
function any_SVE_is_type(
    SV::RetrievalStateVector,
    t::Type{T}
    ) where T<:AbstractStateVectorElement

    for i in eachindex(SV.state_vector_elements)
        if SV.state_vector_elements[i] isa t
            return true
        end
    end

    return false

end


"""
Returns true if SVE is some aerosol-related state vector element. Default is false,
and specific implementations should return true.

# Details
Certain SVEs require dedicated handling in the radiative transfer, in particular during
the computation of atmospheric weighting functions. This function helps us decide if some
SVE should be considered there or not. If an SVE is implmented that *does* require a
dedicated weighting function, users must also provide a function that returns true.
"""
is_aerosol_SVE(SVE::AbstractStateVectorElement) = false


"""
Helper type to assist with selective iteration
over state vector elements that are a subtype
(or type) of `sve_type`.

$(TYPEDFIELDS)

!!! note

    The field `sve_type` holds the type istelf of the desired
    state vector element, and not a variable of the type.

"""
struct StateVectorIterator{T1<:RetrievalStateVector, T2<:AbstractStateVectorElement}

    "Any RetrievalStateVector"
    SV::T1
    "The specific StateVectorElement type to include when iterating"
    sve_type::Type{T2}

end

"""
Implementation of `Base.iterate` to iterate over all
state vector elements in a state vector, but only consider
the type given in `StateVectorIterator`.

Return the tuple `(idx, sve)` where `idx` is the position
of the `StateVectorElement` `sve` in the state vector, if it
is of type `sve_type`.
Very useful, when wanting to loop over state vector elements
of a certain type or type union only. Note that this is
thus not equivalent to iterating with `enumerate`.

$(TYPEDSIGNATURES)


# Usage

To select all state vector elements of type `DispersionPolynomialSVE`,
for example, one can write the following (`my_sv` being an
`AbstractStateVector`).

```
for (idx, sve) in StateVectorIterator(my_sv, DispersionPolynomialSVE)
    println(\"State vector element \$(sve) has index \$(idx).\")
end
```

The above loop is equivalent to writing:

```
for (idx, sve) enumerate(my_sv.state_vector_elements)
    if sve isa DispersionPolynomialSVE
       println(\"State vector element \$(sve) has index \$(idx).\")
    end
end
```

"""
function Base.iterate(s::StateVectorIterator, i::Int=1)

    if i > length(s.SV)
        return nothing
    end

    if s.SV.state_vector_elements[i] isa s.sve_type
        return (i, s.SV.state_vector_elements[i]), i+1
    else
        iterate(s, i+1)
    end
end
