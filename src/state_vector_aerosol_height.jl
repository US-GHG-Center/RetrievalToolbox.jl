"""
Returns the name of this aerosol width state vector element
as a string.

$(SIGNATURES)

"""
function get_name(sve::AerosolHeightSVE)
    if sve.log
        return "AerosolHeightSVE (log) [$(sve.aerosol)]"
    else
        return "AerosolHeightSVE [$(sve.aerosol)]"
    end
end

"""
Pretty printing for AerosolHeightSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::AerosolHeightSVE)

    println(io, "AerosolHeightSVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
Brief pretty printing for AerosolHeightSVE

$(SIGNATURES)
"""
function show(io::IO, sve::AerosolHeightSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end


calculate_jacobian_before_isrf(sve::AerosolHeightSVE) = true
is_aerosol_SVE(sve::AerosolHeightSVE) = true