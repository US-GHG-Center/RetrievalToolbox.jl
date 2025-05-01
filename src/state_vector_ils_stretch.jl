"""
Returns the name of this ILS stretch state vector element
as a string.

$(SIGNATURES)

"""
function get_name(sve::ILSStretchPolynomialSVE)
    return "ILSStretchPolynomialSVE ($(sve.coefficient_order)) " * 
    "[$(sve.swin.window_name)]"
end

"""
Pretty printing for ILSStretchPolynomialSVE type

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::ILSStretchPolynomialSVE)

    println(io, "ILS stretch polynomial SVE")
    println(io, "Window:            $(sve.swin.window_name)")
    println(io, "Coefficient order: $(sve.coefficient_order)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
Brief pretty printing for ILSStretchPolynomialSVE

$(SIGNATURES)
"""
function show(io::IO, sve::ILSStretchPolynomialSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
Returns whether the Jacobian related to ILSStretchPolynomialSVE
types should be calculate before convolution happens.

$(SIGNATURES)
"""
function calculate_jacobian_before_isrf(sve::ILSStretchPolynomialSVE)

    return false

end
