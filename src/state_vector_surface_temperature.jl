"""
Returns the name of this surface temperature state vector element as a string.

$(TYPEDSIGNATURES)

"""
function get_name(sve::SurfaceTemperatureSVE)
    return "SurfaceTemperatureSVE"
end

"""
Pretty printing for `SurfaceTemperatureSVE` types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::SurfaceTemperatureSVE)

    println(io, "Surface temperature SVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
Brief pretty printing for `SurfaceTemperatureSVE`

$(SIGNATURES)
"""
function show(io::IO, sve::SurfaceTemperatureSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `SurfaceTemperatureSVE` types should be calculated
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::SurfaceTemperatureSVE) = true