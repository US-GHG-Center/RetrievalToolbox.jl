"""
Returns the name of this dispersion state vector element
as a string.

$(SIGNATURES)

"""
function get_name(sve::GasLevelScalingFactorSVE)
    return "GasLevelScalingFactorSVE $(sve.gas.gas_name) " *
        "($(sve.start_level), $(sve.end_level))"
end

"""
Pretty printing for DispersionPolynomialSVE types

$(SIGNATURES)
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
Brief pretty printing for DispersionPolynomialSVE

$(SIGNATURES)
"""
function show(io::IO, sve::GasLevelScalingFactorSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

function calculate_jacobian_before_isrf(sve::GasLevelScalingFactorSVE)

    return true

end
