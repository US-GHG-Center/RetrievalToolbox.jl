"""
$(TYPEDSIGNATURES)

Evaluates a polynomial surface at all spectral points indicated by the attached
spectral window. This function is very quick if the surface needs to be evaluated for
all points in the spectral window, as opposed to calling `evalute_surface_at_idx!`
repeatedly.
"""
function evaluate_surface!(
    out::AbstractVector,
    surface::Union{
        LambertianPolynomialSurface,
        LambertianPolynomialKernel,
        RPVPolynomialKernel
    }
)

    swin = surface.swin

    # Quick check
    @assert length(out) == length(swin.ww_grid) "Result array must have the same length " *
        "as the spectral window attached to the surface!"

    poly = ImmutablePolynomial(surface.coefficients)

    # Calculate delta_wavelength = wavelength - ww_reference, where
    # the reference wavelength is the wavelength at the spectral window
    # center.
    delta_wl = swin.ww_grid[:] .- swin.ww_reference

    @turbo for i in eachindex(out)
        out[i] = poly(delta_wl[i])
    end

end


"""
    Evaluates a `LambertianPolynomialSurface` at a given spectral point - this is to be
    used e.g. inside the XRTM module.
"""
function evaluate_surface_at_idx(
    surface::Union{
        LambertianPolynomialSurface,
        LambertianPolynomialKernel,
        RPVPolynomialKernel
        },
    idx::Integer
)

    swin = surface.swin
    poly = ImmutablePolynomial(surface.coefficients)

    delta_wl = swin.ww_grid[idx] - swin.ww_reference

    return poly(delta_wl)

end





"""
$(TYPEDSIGNATURES)

Calculates surface reflectivity due to a Lambertian surface and stores it in-place into
vector `R`.
"""
function calculate_surface_reflectivity!(
    R::AbstractVector,
    surface::LambertianPolynomialSurface,
    scene::EarthScene,
    optical_properties::EarthAtmosphereOpticalProperties
    )
    # Could remove this for a tiny bit of performance gain
    if all(surface.coefficients .== 0)
        @warn "Surface coefficients are all zero!"
    end

    # R will contain the surface albedos
    evaluate_surface!(R, surface)

    # Reflectivity is albedo * µ0 / π
    @views R[:] *= cosd(scene.solar_zenith) / pi

end