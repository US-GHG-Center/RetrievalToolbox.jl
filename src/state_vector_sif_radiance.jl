"""
Returns the name of this temperature offset state vector element as a string.

$(TYPEDSIGNATURES)

"""
function get_name(sve::SIFRadianceSVE)
    return "SIFRadianceSVE"
end

"""
Pretty printing for `SIFRadianceSVE` types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::SIFRadianceSVE)

    println(io, "SIF Radiance SVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
Brief pretty printing for `SIFRadianceSVE`

$(SIGNATURES)
"""
function show(io::IO, sve::SIFRadianceSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `SIFRadianceSVE` types should be calculated
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::SIFRadianceSVE) = true