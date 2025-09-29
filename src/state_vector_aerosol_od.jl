"""
Returns the name of this aerosol OD state vector element as a string.

$(SIGNATURES)

"""
function get_name(sve::AerosolOpticalDepthSVE)
    if sve.log
        return "AerosolOpticalDepthSVE (log) [$(sve.aerosol)]"
    else
        return "AerosolOpticalDepthSVE [$(sve.aerosol)]"
    end
end

"""
Pretty printing for AerosolOpticalDepthSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::AerosolOpticalDepthSVE)

    println(io, "AerosolOpticalDepthSVE [$(sve.aerosol)]")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
Brief pretty printing for AerosolOpticalDepthSVE

$(SIGNATURES)
"""
function show(io::IO, sve::AerosolOpticalDepthSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `AerosolOpticalDepthSVE` types should be calculate
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::AerosolOpticalDepthSVE) = true

"""
$(TYPEDSIGNATURES)

Returns whether this SVE (`sve`) is an aerosol-related SVE. Returns `true`.
"""
is_aerosol_SVE(sve::AerosolOpticalDepthSVE) = true