"""
Returns the name of this dispersion state vector element
as a string.

$(SIGNATURES)

"""
function get_name(sve::DispersionPolynomialSVE)
    return "DispersionPolynomialSVE ($(sve.coefficient_order)) " *
        "[$(sve.dispersion.spectral_window.window_name)]"
end

"""
Pretty printing for DispersionPolynomialSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::DispersionPolynomialSVE)

    println(io, "Dispersion polynomial SVE")
    println(io, "Coefficient order: $(sve.coefficient_order)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
Brief pretty printing for DispersionPolynomialSVE

$(SIGNATURES)
"""
function show(io::IO, sve::DispersionPolynomialSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
Returns whether the Jacobian related to DispersionPolynomialSVE
types should be calculate before convolution happens.

$(SIGNATURES)
"""
function calculate_jacobian_before_isrf(sve::DispersionPolynomialSVE)

    return false

end
