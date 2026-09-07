"""
Returns the name of this `ZeroLevelOffsetSVE` state vector element as a string.

$(TYPEDSIGNATURES)

"""
function get_name(sve::ZeroLevelOffsetPolynomialSVE)
    return "ZeroLevelOffsetPolynomialSVE ($(sve.coefficient_order)) " *
        "[$(sve.swin.window_name)]"
end

"""
Pretty printing for ZeroLevelOffsetSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::ZeroLevelOffsetPolynomialSVE)

    println(io, "Zero level offset polynomial SVE")
    println(io, "Coefficient order: $(sve.coefficient_order)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Current value:     $(sve.iterations[end])")

end

"""
Brief pretty printing for SurfacePressureSVE

$(SIGNATURES)
"""
function show(io::IO, sve::ZeroLevelOffsetPolynomialSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `ZeroLevelOffsetPolynomialSVE` types should be
calculated before convolution happens. Returns `false` since the partial derivative
calculation does not need the resulting radiance itself.
"""
calculate_jacobian_before_isrf(sve::ZeroLevelOffsetPolynomialSVE) = false