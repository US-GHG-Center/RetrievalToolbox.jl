"""
$(TYPEDSIGNATURES)

Returns the name of this surface albedo state vector element as a string.
"""
function get_name(sve::SurfaceAlbedoPolynomialSVE)
    return "SurfaceAlbedoPolynomialSVE ($(sve.coefficient_order)) " *
        "[$(sve.swin.window_name)]"

end


"""
$(TYPEDSIGNATURES)

Pretty printing for `SurfaceAlbedoPolynomialSVE` types
"""
function show(io::IO, ::MIME"text/plain", sve::SurfaceAlbedoPolynomialSVE)

    println(io, "Surface albedo polynomial SVE")
    println(io, "Window:            $(sve.swin.window_name)")
    println(io, "Coefficient order: $(sve.coefficient_order)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for `SurfaceAlbedoPolynomialSVE`
"""
function show(io::IO, sve::SurfaceAlbedoPolynomialSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `SurfaceAlbedoPolynomialSVE` types should be
calculated before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::SurfaceAlbedoPolynomialSVE) = true