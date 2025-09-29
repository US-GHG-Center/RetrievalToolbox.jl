"""
$(TYPEDSIGNATURES)

Returns the name of this surface albedo state vector element as a string.
"""
function get_name(sve::BRDFPolynomialSVE)
    BRDF_type_short = get_short_name(sve.BRDF_type)
    return "BRDFPolynomialSVE $(BRDF_type_short) ($(sve.coefficient_order)) " *
        "[$(sve.swin.window_name)]"

end


"""
$(TYPEDSIGNATURES)

Pretty printing for `BRDFPolynomialSVE` types
"""
function show(io::IO, ::MIME"text/plain", sve::BRDFPolynomialSVE)

    println(io, "Surface BRDF polynomial SVE")
    println(io, "Window:            $(sve.swin.window_name)")
    println(io, "Coefficient order: $(sve.coefficient_order)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for `BRDFPolynomialSVE`
"""
function show(io::IO, sve::BRDFPolynomialSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `BRDFPolynomialSVE` types should be calculate
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::BRDFPolynomialSVE) = true
