"""
$(TYPEDSIGNATURES)

Get the isotropic SIF radiance for a requested wavelength. Must supply a wavelength or
wavenumber *with* units. NOTE: this is a rather slow function. For better performance,
use `get_SIF_radiance!` when you need to evaluate SIF radiance for a whole vector of
spectral points.
"""
function get_SIF_radiance(
    sif::SIFRadiance,
    ww::Union{Unitful.Length,Unitful.Wavenumber}
    )


    # Convert ww into the same units as the SIF radiance
    if unit(ww) isa Unitful.WavenumberUnits
        # Must convert to wavelength, since the SIF waveform is natively stored in
        # wavelength domain.
        wl = 1 / ww |> sif.ww_unit
    else
        wl = ww |> sif.ww_unit
    end

    if (wl < sif.ww_waveform[1]) || (wl > sif.ww_waveform[end])
        return 0.0 # zero is zero regardless of units
    end

    # wl is now the user-requested spectral point in units of sif.ww_unit

    #=
        Return the SIF radiance at wavelength wl, and make the appropriate scaling:
        sif.radiance_at_reference tells us what the SIF radiance at the reference spectral
        point is, and we have to scale the waveform
    =#

    return sif.radiance_at_reference * sif.itp(wl) /
        sif.itp(sif.ww_reference * sif.ww_unit)

end


"""
$(TYPEDSIGNATURES)

"""
function get_SIF_radiance!(
    result::Vector{<:Number},
    sif::SIFRadiance,
    ww_grid::Vector{<:Number},
    ww_unit::Unitful.WavenumberUnits
)
    @assert length(result) == length(ww_grid)

    # Convert ww into the same units as the SIF radiance
    unit_fac_ww = 1 / ww_unit |> sif.ww_unit

    # SIF scale factor
    SIF_scale = sif.radiance_at_reference * sif.itp(sif.ww_reference * sif.ww_unit)

    for i in eachindex(ww_grid)

        wl = unit_fac_ww / ww_grid[i]# This in wavelengths now
        result[i] = SIF_scale * sif.itp(wl)

    end

end

"""
$(TYPEDSIGNATURES)

"""
function get_SIF_radiance!(
    result::Vector{<:Number},
    sif::SIFRadiance,
    ww_grid::Vector{<:Number},
    ww_unit::Unitful.LengthUnits
)
    @assert length(result) == length(ww_grid)

    # Convert ww into the same units as the SIF radiance
    unit_fac_ww = ww_unit |> sif.ww_unit

    # SIF scale factor
    SIF_scale = sif.radiance_at_reference * sif.itp(sif.ww_reference * sif.ww_unit)

    for i in eachindex(ww_grid)

        wl = unit_fac_ww * ww_grid[i] # This in wavelengths now
        result[i] = SIF_scale * sif.itp(wl)

    end

end