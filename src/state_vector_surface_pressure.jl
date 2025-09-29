"""
Returns the name of this surface pressure vector element as a string.

$(TYPEDSIGNATURES)

"""
function get_name(sve::SurfacePressureSVE)
    return "SurfacePressureSVE"
end

"""
Pretty printing for `SurfacePressureSVE` types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::SurfacePressureSVE)

    println(io, "Surface pressure SVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
Brief pretty printing for `SurfacePressureSVE`

$(SIGNATURES)
"""
function show(io::IO, sve::SurfacePressureSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `SurfacePressureSVE` types should be calculated
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::SurfacePressureSVE) = true