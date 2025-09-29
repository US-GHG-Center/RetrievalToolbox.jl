"""
$(TYPEDSIGNATURES)

Returns the name of this aerosol width state vector element as a string.
"""
function get_name(sve::AerosolHeightSVE)
    if sve.log
        return "AerosolHeightSVE (log) [$(sve.aerosol)]"
    else
        return "AerosolHeightSVE [$(sve.aerosol)]"
    end
end

"""
$(TYPEDSIGNATURES)

Pretty printing for `AerosolHeightSVE` types
"""
function show(io::IO, ::MIME"text/plain", sve::AerosolHeightSVE)

    println(io, "AerosolHeightSVE")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Last iteration:    $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for `AerosolHeightSVE`
"""
function show(io::IO, sve::AerosolHeightSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end


"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `AerosolHeightSVE` types should be calculate
before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::AerosolHeightSVE) = true

"""
$(TYPEDSIGNATURES)

Returns whether this SVE (`sve`) is an aerosol-related SVE. Returns `true`.
"""
is_aerosol_SVE(sve::AerosolHeightSVE) = true