"""
$(TYPEDSIGNATURES)

Returns the name of this dispersion state vector element as a string.
"""
function get_name(sve::GasLevelScalingFactorSVE)
    return "GasLevelScalingFactorSVE $(sve.gas.gas_name) " *
        "($(sve.start_level), $(sve.end_level))"
end

"""
$(TYPEDSIGNATURES)

Pretty printing for `DispersionPolynomialSVE` types
"""
function show(io::IO, ::MIME"text/plain", sve::GasLevelScalingFactorSVE)

    println(io, "GasLevelScalingFactorSVE")
    println(io, "Start level:       $(sve.start_level)")
    println(io, "End level:         $(sve.end_level)")
    println(io, "Corresponding gas: $(sve.gas.gas_name)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for `GasLevelScalingFactorSVE`.
"""
function show(io::IO, sve::GasLevelScalingFactorSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `GasLevelScalingFactorSVE` types should be
calculated before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::GasLevelScalingFactorSVE) = true