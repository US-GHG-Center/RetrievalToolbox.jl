"""
$(TYPEDSIGNATURES)

Converts a symbol to its corresponding one that carries the unit information. E.g.
`:irradiance` returns `:irradiance_unit`, or `:temperature_levels` returns
`:temperature_unit`
"""
function _field_unit_conversion(x::Symbol)

    return Symbol(replace(String(x), "_layers" => "", "_levels" => "") * "_unit")

end

"""
$(TYPEDSIGNATURES)

Sets a RetrievalToolbox type field and applies the correct unit conversion

# Details
"""
function ingest!(
    obj,
    field::Symbol,
    val::Number
)

    # Create the symbol that references the field unit
    # (e.g. :temperature_unit to :temperature)
    obj_unit_symbol = _field_unit_conversion(field)

    # Check if the object has the field
    if !hasproperty(obj, obj_unit_symbol)
        @error "Object $(obj) does not have $(field)!"
    end

    # Grab the object
    obj_field = getproperty(obj, field)
    # Grab the unit
    obj_unit = getproperty(obj, obj_unit_symbol)

    # If the object field is a number, we can simply do this:
    if obj_field isa Number
        # Set the field with the correct units - note that we have
        # to strip the units to assign the "bare" quantity.
        setfield!(obj, field, val |> obj_unit |> ustrip)
        return nothing
    # If the object field is some array that we want to set to the same value,
    # do the following:
    elseif obj_field isa AbstractArray
        # Pre-compute what we want to insert
        rhs = val |> obj_unit .|> ustrip
        # Loop to insert values
        for i in eachindex(obj_field)
            obj_field[i] = rhs
        end
    end

end

"""
$(TYPEDSIGNATURES)

Sets a RetrievalToolbox type field and applies the correct unit conversion. This is the specialized
function to deal with array-type values. This should be mostlya allocation-free.

# Details
"""
function ingest!(
    obj,
    field::Symbol,
    val::AbstractArray
)

    # Get the reference to the array from the object
    obj_field = getproperty(obj, field)
    # Create the symbol that references the field unit
    # (e.g. :temperature_unit to :temperature)
    obj_unit_symbol = _field_unit_conversion(field)
    # Grab the unit
    obj_unit = getproperty(obj, obj_unit_symbol)

    # Calculate unit conversion factor, so this
    # will be a number that the @turbo macro below can interpret
    uc = one(eltype(val)) * unit(val[1]) |> obj_unit |> ustrip

    # Get the unit-stripped view to the "val" array
    ustrip_val = ustrip.(val)

    # Copy over, accounting for correct conversion factor
    for i in eachindex(ustrip_val)
        obj_field[i] = ustrip_val[i] * uc
    end
end

"""
$(TYPEDSIGNATURES)

For an array `x`, returns the difference between the largest and smallest element
`maximum(x) - minimum(x)`.
"""
function maxmin(x::AbstractArray)
    return maximum(x) - minimum(x)
end

"""
$(TYPEDSIGNATURES)

For an array `x`, returns `true` if any element is not finite (as per `isfinite`), or
`false` if all elements are finite.
"""
function check_for_not_finite(x::AbstractArray)

    for i in eachindex(x)
        if !isfinite(x[i])
            return true
        end
    end

    return false

end

function avx_sum_along_columns_between!(y, x, idx1, idx2)

    @turbo for j in idx1:idx2
        for i in axes(x, 1)

            y[i] += x[i,j]

        end
    end

end


function avx_add_along_columns!(y, x)
    @turbo for j in axes(x, 2)
        for i in axes(x, 1)

            y[i] += x[i,j]

        end
    end
end

function avx_sum_along_columns!(y, x)
    @views y[:] .= 0.0
    @turbo for j in axes(x, 2)
        for i in axes(x, 1)

            y[i] += x[i,j]

        end
    end
end

function avx_sum(x)

    # Initialize this value to be of the input
    total = zero(eltype(x))

    @turbo for i in eachindex(x)
        total += x[i]
    end

    return total

end

function avx_dot(x::AbstractVector,
                 y::AbstractVector)

    @assert length(x) == length(y) "Both vectors must be the same length"

    # Infer the output variable type from "x"
    total = zero(eltype(x))

    @turbo for i in eachindex(x)
        total += x[i] * y[i]
    end

    return total

end


function vector_contained_in_vector(vec1::AbstractVector,
                                    vec2::AbstractVector)

    len1 = length(vec1)
    len2 = length(vec2)

    if len2 > len1
        return -1
    end

    for i in 1:len1 - len2 + 1
        if @views vec1[i:i+len2-1] == vec2[:]
            return i
        end
    end

    # Nothing found .. return -1
    return -1

end

"""
$(TYPEDSIGNATURES)

For a given level object, create a new vector that contains the mid-layer vector. This
function assumes that the layer-value is evaluated at the center point between the two
adjacent levels.
"""
function levels_to_layers(levels::AbstractVector; logspace=false)

    # Strip units here, taking log is a little counterintuitive
    # for the Uniful package.
    N_levels = length(levels)
    layers = zeros(eltype(levels), N_levels - 1)

    for i in 1:N_levels - 1
        layers[i] = 0.5 * (levels[i] + levels[i + 1])
    end

    return layers

end

"""
$(TYPEDSIGNATURES)

For a given `level` vector, calculate the mid-layer values and store into `layers`. This
function assumes that the layer-value is evaluated at the center point between the two
adjacent levels.
"""
function levels_to_layers!(
    layers::AbstractVector,
    levels::AbstractVector;
    logspace=false
    )

    # Strip units here, taking log is a little counterintuitive
    # for the Uniful package.
    N_levels = length(levels)

    for i in 1:N_levels - 1
        layers[i] = 0.5 * (levels[i] + levels[i + 1])
    end

end


function atmospheric_profile_interpolator_linear(
    plevels_hires::AbstractVector,
    atm_hires::AbstractVector,
    plevels_out::AbstractVector;
    logspace=false,
    cutoff_threshold=0.05
    )

    # Strip units here, taking log is a little counterintuitive
    # for the Uniful package, but convert to output units first
    _unit_plevels_out = unit(plevels_out[1])
    _plevels_hires = ustrip(uconvert.(_unit_plevels_out, plevels_hires))
    _plevels_out = ustrip(plevels_out)

    if logspace
        _plevels_hires = log.(_plevels_hires)
        _plevels_out = log.(_plevels_out)
    end

    atm_out = LinearInterpolation(
        _plevels_hires, atm_hires, extrapolation_bc=Linear())(_plevels_out)

    return atm_out

end


function pwl_value_1d(xd, yd, xi)

    last_d = 1

    ni = length(xi)
    nd = length(xd)
    yi = zeros(eltype(yd), ni)

    idx_inside_left = 1
    for i in eachindex(xi)
        if xi[i] <= xd[1]
            yi[i] = yd[1]
            idx_inside_left = i + 1
        else
            break
        end
    end

    for i in idx_inside_left:idx_inside_right
        for d in last_d:nd-1
            if (xi[i] >= xd[d]) & (xi[i] <= xd[d+1])
                last_d = d
                @inbounds fac = (xi[i] - xd[d]) / (xd[d+1] - xd[d])
                @inbounds yi[i] = (1.0 - fac) * yd[d] + fac * yd[d+1]
                break
            end
        end
    end

    return yi

end

function pwl_value_1d!(xd, yd, xi, yi)

    last_d = 1

    nd = length(xd)

    idx_inside_left = 1
    for i in eachindex(xi)
        if xi[i] <= xd[1]
            yi[i] = yd[1]
            idx_inside_left = i + 1
        else
            break
        end
    end


    idx_inside_right = length(xi)
    for i in length(xi):-1:1
        if xi[i] >= xd[end]
            yi[i] = yd[end]
            idx_inside_right = i - 1
        else
            break
        end
    end


    for i in idx_inside_left:idx_inside_right
        for d in last_d:nd-1
            if (xi[i] >= xd[d]) & (xi[i] <= xd[d+1])
                last_d = d
                @inbounds fac = (xi[i] - xd[d]) / (xd[d+1] - xd[d])
                @inbounds yi[i] = (1.0 - fac) * yd[d] + fac * yd[d+1]
                break
            end
        end
    end

end

"""
a * b

"""
function pwl_value_1d_axb!(xd, ad, bd, xi, yi)

    last_d = 1

    nd = length(xd)

    idx_inside_left = 1
    for i in eachindex(xi)
        if xi[i] <= xd[1]
            yi[i] = ad[1] * bd[1]
            idx_inside_left = i + 1
        else
            break
        end
    end


    idx_inside_right = length(xi)
    for i in length(xi):-1:1
        if xi[i] >= xd[end]
            yi[i] = ad[end] * bd[end]
            idx_inside_right = i - 1
        else
            break
        end
    end


    for i in idx_inside_left:idx_inside_right
        for d in last_d:nd-1
            if (xi[i] >= xd[d]) & (xi[i] <= xd[d+1])
                last_d = d
                @inbounds fac = (xi[i] - xd[d]) / (xd[d+1] - xd[d])
                @inbounds yi[i] = (1.0 - fac) * ad[d] * bd[d] + fac * ad[d+1] * bd[d+1]
                break
            end
        end
    end

end

"""
$(TYPEDSIGNATURES)

Calculates the trapezoidal-rule integral of f(x) where `x` are discrete x_i, and `y` are
the corresponding f(x_i).
"""
function _trapz(x, y)

    result = zero(eltype(y))

    for i in 2:length(y)
        result += 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    end

    return result

end

"""
$(TYPEDSIGNATURES)

Returns `true` if `a` is found within `x`, does not allocate.

# Examples
```jldoctest
julia> findany(["a", "b", "c", "d"], "b")
true
```
```jldoctest
julia> findany([1,2,3,4,5], 10)
false
```
"""
function findany(x::AbstractVector, a)
    @inbounds for i in eachindex(x)
        x[i] == a && return true
    end
    return false
end


"""
$(TYPEDSIGNATURES)

Returns `true` if any element within `x` is of type `T`, does not allocate.
"""
function findanytype(x::AbstractVector, T)
    @inbounds for i in eachindex(x)
        x[i] isa T && return true
    end
    return false
end

"""
$(TYPEDSIGNATURES)

Calculates specific humidity from H2O VMR
"""
function H2O_VMR_to_specific_humidity(
    h2o::Union{Number, Unitful.DimensionlessQuantity}
)
    # H2O can have units, but we strip them here
    x = h2o |> NoUnits

    # Specific humidity calculation, returned in
    # [kg/kg] or [1]
    q = h2o / ((1 - h2o) * MM_AIR_TO_H2O + h2o)

    return q
end

"""
$(TYPEDSIGNATURES)

Calculates H2O VMR from specific humidity
"""
function specific_humidity_to_H2O_VMR(
    sh::Union{Number, Unitful.DimensionlessQuantity}
)
    # H2O can have units, but we strip them here
    # (this automatically converts)
    q = sh |> NoUnits

    # Specific humidity calculation, returned in
    # [kg/kg] or [1]
    h2o  = q / ((1 - q) * MM_H2O_TO_AIR + q)

    return h2o
end

"""
$(TYPEDSIGNATURES)

Calculates the standard deviation of a Gaussian for a given full width at half the maximum
value.
"""
function FWHM_to_sigma(FWHM::Number)

    return FWHM / 2 / (sqrt(2 * log(2)))

end