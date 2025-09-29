"""
$(TYPEDSIGNATURES)

Returns the name of this aerosol width state vector element as a string.
"""
function get_name(sve::AerosolWidthSVE)
    if sve.log
        return "AerosolWidthSVE (log) [$(sve.aerosol)]"
    else
        return "AerosolWidthSVE [$(sve.aerosol)]"
    end
end

"""
$(TYPEDSIGNATURES)

Pretty printing for AerosolWidthSVE types
"""
function show(io::IO, ::MIME"text/plain", sve::AerosolWidthSVE)

    println(io, "AerosolWidthSVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for AerosolWidthSVE
"""
function show(io::IO, sve::AerosolWidthSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `AerosolWidthSVE` types should be calculate
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::AerosolWidthSVE) = true

"""
$(TYPEDSIGNATURES)

Returns whether this SVE (`sve`) is an aerosol-related SVE. Returns `true`.
"""
is_aerosol_SVE(sve::AerosolWidthSVE) = true