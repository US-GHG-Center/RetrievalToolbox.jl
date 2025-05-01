"""
Returns the name of this temperature offset state vector element
as a string.

$(SIGNATURES)

"""
function get_name(sve::TemperatureOffsetSVE)
    return "TemperatureOffsetSVE"
end

"""
Pretty printing for TemperatureOffsetSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::TemperatureOffsetSVE)

    println(io, "Temperature Offset SVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
Brief pretty printing for TemperatureOffsetSVE

$(SIGNATURES)
"""
function show(io::IO, sve::TemperatureOffsetSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end


function calculate_jacobian_before_isrf(sve::TemperatureOffsetSVE)

    return true

end
