#=

    General Notes

    Beer-Lambert type RT calculations will only ever produce radiance intensities, so
    everything below is hardcoded to only modify the intensity component (.I). However,
    there is no explicit restriction to `ScalarRadaiance` types in these functions below,
    so users could provide `VectorRadiance` objects in the corresponding RT structures,
    but the other componentes (Q,U) will not be calculated or modified.

=#


"""
$(TYPEDSIGNATURES)

Pretty printing for Beer-Lambert RT
"""
function show(io::IO, ::MIME"text/plain", rt::BeerLambertRTMethod)
    # Just print this for now
    println(io, "Beer-Lambert RT")

end

"""
$(TYPEDSIGNATURES)

Pretty printing for Beer-Lambert RT
"""
function show(io::IO, rt::BeerLambertRTMethod)
    # Just print this for now
    println(io, "Beer-Lambert RT")

end


"""
$(TYPEDSIGNATURES)

Calculates radiances and Jacobians for a `BeerLambertRTMethod` object. This
further dispatches to the correct function for the specific observer type.
"""
function calculate_radiances_and_jacobians!(rt::BeerLambertRTMethod)

    # Make explicit dispatch to function, depending on observer mode
    calculate_radiances_and_jacobians!(rt, rt.scene.observer)

end

"""
$(TYPEDSIGNATURES)

Calculates the radiances and Jacobians of a `BeerLambertRTMethod` object,
for a `SatelliteObserver` type.
"""
function calculate_radiances_and_jacobians!(
    rt::BeerLambertRTMethod,
    observer::SatelliteObserver
    )

    # Rebind to surface for faster access
    surface = rt.scene.surfaces[rt.optical_properties.spectral_window]

    # Calculate surface reflectivity
    # inside temporary array.
    R = rt.optical_properties.tmp_Nhi1
    @views R[:] .= 0.0

    # This shall dispatch to the specific
    # surface type supplied in the scene.
    calculate_surface_reflectivity!(
        R,
        surface,
        rt.scene,
        rt.optical_properties
    )

    # Total column optical depth
    total_column_od = rt.optical_properties.tmp_Nhi2
    @views total_column_od[:] .= 0.0

    # Zero out hires radiance container
    @views rt.hires_radiance.I[:] .= 0.0

    # In-place summation, equivalent to
    # @views total_column_od[:] = sum(rt.optical_properties.total_tau, dims=2)[:,1]
    # but non-allocating, and overwriting the contents of total_column_od
    avx_sum_along_columns!(total_column_od, rt.optical_properties.total_tau)
    rt.hires_radiance.I[:] = total_column_od[:]

    # At this point, rt.hires_radiance would be calculated in units
    # rt.scene.solar_model.unit / u"sr"
    # Any Jacobian of a state vector element "sve" has units of
    # rt.scene.solar_model.unit / u"sr" / get_unit(sve)

    # Top of the atmosphere radiance is simply
    # the incoming beam multiplied by the surface reflectance
    # and modified by the optical extinction of the atmosphere
    # for both incoming and outgoing paths.
    @turbo for i in eachindex(rt.hires_radiance.I)

        rt.hires_radiance.I[i] = (
            rt.hires_solar.I[i] * R[i] * exp(
                -total_column_od[i] / cosd(rt.scene.solar_zenith)
                -total_column_od[i] / cosd(rt.scene.observer.viewing_zenith)
            )
        )

        # Due to the unit factor, however, we scale up the radiance units
        # to be of the same unit as demanded by `rt.unit`.

    end

    #= NOTE
    cosd(10.0) == cosd(10.0u"°")
    cosd/sind calculate sin/cos for a degree-valued argument, and if the argument
    is a Unitful quantity in degrees, it will act accordingly.
    HOWEVER, one must not supply e.g. cosd(10.0u"rad"), as NO FURTHER unit conversion
    will take place, so
    cosd(10.0) == cosd(10.0u"°") != cosd(10.0u"rad")
    =#


    #=
        NOTE
        Units in rt.hires_radiance and rt.hires_jacobians might not be the
        same as in the RT buffer object! This needs to be pointed out clearly
        in the documentation: the solar model governs the units used in the
        RT object (not the RT buffer though!)
    =#


    # Calculate jacobians
    for (i, sve) in enumerate(rt.state_vector.state_vector_elements)
        if calculate_jacobian_before_isrf(sve)
            calculate_rt_jacobian!(rt.hires_jacobians[sve], rt, sve)
        end
    end


end

function calculate_radiances_and_jacobians!(
    rt::BeerLambertRTMethod,
    observer::UplookingGroundObserver
    )


    # Calculate sphericity factors
    ch = create_sphericity_factors(rt.optical_properties, rt.scene)
    #sza_per_lay = create_refracted_sza(rt.optical_properties, rt.scene)

    # Total column optical depth
    total_column_od = rt.optical_properties.tmp_Nhi2
    @views total_column_od[:] .= 0.0

    @views rt.hires_radiance.I[:] .= 0.0

     # In-place summation, equivalent to
    # @views total_column_od[:] = sum(rt.optical_properties.total_tau, dims=2)[:,1]
    # but non-allocating, and overwriting the contents of total_column_od
    #avx_sum_along_columns!(total_column_od, rt.optical_properties.total_tau)


    @turbo for l in 1:rt.scene.atmosphere.N_layer
        for i in eachindex(rt.hires_radiance.I)
            #total_column_od[i] += rt.optical_properties.total_tau[i,l] / cosd(rt.scene.solar_zenith)
            total_column_od[i] += rt.optical_properties.total_tau[i,l] * ch[l]
            #total_column_od[i] += rt.optical_properties.total_tau[i,l] / cosd(sza_per_lay[l])
        end
    end


    # At this point, rt.hires_radiance would be calculated in units
    # rt.scene.solar_model.unit / u"sr"
    # Any Jacobian of a state vector element "sve" has units of
    # rt.scene.solar_model.unit / u"sr" / get_unit(sve)

    # Same as for a SatelliteObserver, but for the UplookingObserver, we
    # drop the outgoing path as well as the surface reflectance.
    @turbo for i in eachindex(rt.hires_radiance.I)
        rt.hires_radiance.I[i] = rt.hires_solar.I[i] * exp(-total_column_od[i]
            #/ cosd(rt.scene.solar_zenith)
            )
    end


   # Calculate jacobians
    for (i, sve) in enumerate(rt.state_vector.state_vector_elements)
        if calculate_jacobian_before_isrf(sve)
            calculate_rt_jacobian!(rt.hires_jacobians[sve], rt, sve)
        end
    end


end

"""
$(TYPEDSIGNATURES)

Default behavior for `BeerLambertRTMethod` jacobians
"""
function calculate_rt_jacobian!(
    rt::BeerLambertRTMethod,
    sve::AbstractStateVectorElement
)
    @debug "[RT] No Beer-Lambert RT Jacobian calculated for $(sve)"
    return nothing
end


"""
$(TYPEDSIGNATURES)

Calculates the `SurfaceAlbedoPolynomialSVE` Jacobian for a `BeerLambertRTMethod` RT
object. This further dispatches to the correct function for the specific observer type.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::SurfaceAlbedoPolynomialSVE
)

    # Dispatch to specific observer type
    calculate_rt_jacobian!(jac, rt, sve, rt.scene.observer)

end


"""
$(TYPEDSIGNATURES)

Calculates the `SurfaceAlbedoPolynomialSVE` Jacobian for a `BeerLambertRTMethod` RT
object and a `SatelliteObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::SurfaceAlbedoPolynomialSVE,
    observer::SatelliteObserver
    )

    # Rebind to surface for more convenient access
    surface = rt.scene.surfaces[rt.optical_properties.spectral_window]
    swin = rt.optical_properties.spectral_window

    # Simply return if this RT method object is not the one
    # the SVE is referring to!
    # E.g., a surface albedo SVE that is attached to Band X
    # should report 0.0 for Band Y
    if !(sve.swin === rt.optical_properties.spectral_window)
        @views jac.I[:] .= 0.0
        return nothing
    end

    # Calculate dI / dAlbedo_coefficient
    # dI / dAlbedo_coefficient = dI / dAlbedo * dAlbedo / dAlbedo_coefficient
    #                          = I / Albedo * (wavelength - ww_reference)^order

    # Keep in mind that R = albedo * mu0 / pi, so that
    # mu0 / pi factor needs to come back later.

    # Zero out
    R = rt.optical_properties.tmp_Nhi1
    @views R[:] .= 0.0

    # Calculate R
    calculate_surface_reflectivity!(
        R,
        surface,
        rt.scene,
        rt.optical_properties
    )

    @views @. jac.I[:] = rt.hires_radiance.I[:] / R[:] *
        (cosd(rt.scene.solar_zenith) / pi)

    if sve.coefficient_order > 0

        delta_wl = swin.ww_grid[:] .- swin.ww_reference

        @turbo for i in eachindex(jac.I)
            jac.I[i] *= (delta_wl[i] ^ (sve.coefficient_order))
        end

    end

    # zero out tmps
    @views rt.optical_properties.tmp_Nhi1[:] .= 0.0
    @views rt.optical_properties.tmp_Nhi2[:] .= 0.0

end

"""
$(TYPEDSIGNATURES)

Calculates the `SurfacePressureSVE` Jacobian for a `BeerLambertRTMethod` RT object. This
further dispatches to the correct function for the specific observer type.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::SurfacePressureSVE,
    )

    # Dispatch to specific observer type
    calculate_rt_jacobian!(jac, rt, sve, rt.scene.observer)

end

"""
$(TYPEDSIGNATURES)

Calculates the `SurfacePressureSVE` Jacobian for a `BeerLambertRTMethod` RT object and a
`SatelliteObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::SurfacePressureSVE,
    observer::SatelliteObserver
    )

    # What gases do we have with dTau/dpsurf?
    N_hires = rt.optical_properties.spectral_window.N_hires
    gases = keys(rt.optical_properties.gas_derivatives)

    # Unit factor needed to get from SVE unit to Pa
    # E.g., if the surface pressure SVE is in hPa,
    # the calculated Jacobian must be adjusted by a factor
    # 1 hPa / 1 Pa = 100.0
    unit_factor = ustrip(rt.scene.atmosphere.pressure_unit, 1.0 * get_unit(sve))

    # Take a buffer array, and zero out
    dtau_dpsurf = rt.optical_properties.tmp_Nhi1
    @views @. dtau_dpsurf[:] = 0

    # Sum up ∂Tau/∂psurf for all gases
    for gas in gases
        @views @. dtau_dpsurf[:] +=
            rt.optical_properties.gas_derivatives[gas]["dTau_dpsurf"][:]
    end

    # Equation:
    # dI / dpsurf = dI / dTau * dTau / dpsurf
    #             = - I (1/mu_0 + 1/mu) * dTau/dpsurf
    @turbo for i in eachindex(jac.I)

        jac.I[i] = -rt.hires_radiance.I[i] * (
            1 / cosd(rt.scene.solar_zenith) +
            1 / cosd(rt.scene.observer.viewing_zenith)
        ) * dtau_dpsurf[i] * unit_factor
    end

end

"""
$(TYPEDSIGNATURES)

Calculates the `SurfacePressureSVE` Jacobian for a `BeerLambertRTMethod` RT object and a
`UplookingGroundObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::SurfacePressureSVE,
    observer::UplookingGroundObserver
    )

    # What gases do we have with dTau/dpsurf?
    N_hires = rt.optical_properties.spectral_window.N_hires
    gases = keys(rt.optical_properties.gas_derivatives)

    # Unit factor needed to get from SVE unit to Pa
    # E.g., if the surface pressure SVE is in hPa,
    # the calculated Jacobian must be adjusted by a factor
    # 1 hPa / 1 Pa = 100.0
    unit_factor = ustrip(rt.scene.atmosphere.pressure_unit, 1.0 * get_unit(sve))

    # Take a buffer array, and zero out
    dtau_dpsurf = rt.optical_properties.tmp_Nhi1
    @views @. dtau_dpsurf[:] = 0

    # Sum up ∂Tau/∂psurf for all gases
    for gas in gases
        @views @. dtau_dpsurf[:] +=
            rt.optical_properties.gas_derivatives[gas]["dTau_dpsurf"][:]
    end

    # Equation:
    # dI / dpsurf = dI / dTau * dTau / dpsurf
    #             = - I (1/mu_0) * dTau/dpsurf
    @turbo for i in eachindex(jac.I)

        jac.I[i] = -rt.hires_radiance.I[i] * (
            1 / cosd(rt.scene.solar_zenith)
        ) * dtau_dpsurf[i] * unit_factor
    end

end


"""
$(TYPEDSIGNATURES)

Calculates the `GasVMRProfileSVE` Jacobian for a `BeerLambertRTMethod` RT object. This
further dispatches to the correct function for the specific observer type.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::GasVMRProfileSVE,
    )

    @views jac.I[:] .= 0

    # If this gas is not part of this window, ignore
    if !haskey(rt.optical_properties.gas_tau, sve.gas)
        return
    end

    # Dispatch to specific observer type
    calculate_rt_jacobian!(jac, rt, sve, rt.scene.observer)

end


function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::GasVMRProfileSVE,
    observer::SatelliteObserver
    )

    # Check if the gas from the SVE is actually present
    if !haskey(rt.optical_properties.gas_derivatives, sve.gas)
        @error "Gas $(sve.gas) not present in optical properties of $(rt)"
        return nothing
    end

    # Use tmp array for this calculation
    dTau_dlevel = rt.optical_properties.tmp_Nhi1
    # Zero-out, just to be sure..
    @views dTau_dlevel[:] .= 0.0
    # How many levels does our retrieval grid have?
    N_level = rt.scene.atmosphere.N_level

    # Bind ∂τ/∂VMR for ease of access
    dTau_dVMR = rt.optical_properties.gas_derivatives[sve.gas]["dTau_dVMR"]

    # Depending on the level, we construct the ∂τ_layer / ∂VMR_level differently
    if sve.level == 1
        @turbo for i in eachindex(dTau_dlevel)
            dTau_dlevel[i] = dTau_dVMR[i,sve.level,1]
        end
    elseif sve.level == N_level
        @turbo for i in eachindex(dTau_dlevel)
            dTau_dlevel[i] = dTau_dVMR[i,sve.level-1,2]
        end
    else
        @turbo for i in eachindex(dTau_dlevel)
            dTau_dlevel[i] = dTau_dVMR[i,sve.level-1,2] + dTau_dVMR[i,sve.level,1]
        end
    end

    # In dTau_dVMR, the dVMR is per unitless quantity,
    # so we must introduce a conversion factor to account
    # for the difference in the state vector unit (e.g. ppm or ppb)
    # to this "per VMR in units of one"
    unit_fac = ustrip(sve.unit, 1.0)

    @turbo for i in eachindex(jac.I)

        jac.I[i] = -rt.hires_radiance.I[i] * (
            (1.0 / cosd(rt.scene.solar_zenith)) +
                (1.0 / cosd(rt.scene.observer.viewing_zenith))
        ) * dTau_dlevel[i] / unit_fac

        # Note the final division by the unit factor which makes sure that the final
        # quantity is indeed in units of SVE_unit / radiance_unit.

    end

    # Re-zero!
    @views dTau_dlevel[:] .= 0.0

end


"""
$(TYPEDSIGNATURES)


Calculates the `GasLevelScalingFactorSVE` Jacobian for a `BeerLambertRTMethod` RT object. This
further dispatches to the correct function for the specific observer type.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::GasLevelScalingFactorSVE,
    )

    @views jac.I[:] .= 0

    # If this gas is not part of this window, ignore
    if !haskey(rt.optical_properties.gas_tau, sve.gas)
        return
    end

    # Dispatch to specific observer type
    calculate_rt_jacobian!(jac, rt, sve, rt.scene.observer)

end

"""
$(TYPEDSIGNATURES)

Calculates the `GasLevelScalingFactorSVE` Jacobian for a `BeerLambertRTMethod` RT object
and a `SatelliteObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::GasLevelScalingFactorSVE,
    observer::SatelliteObserver
    )

    # If this gas is not part of this window, ignore
    if !haskey(rt.optical_properties.gas_tau, sve.gas)
        @debug "[RT] Skipping RT Jacobian for $(sve.gas)"
        return
    end

    # current value of scale factor, and apply units
    scale_factor = sve.iterations[end] * sve.unit |> NoUnits

    # Layer summation start and end
    idx1 = sve.start_level
    idx2 = sve.end_level - 1

    tau_subcolumn = rt.optical_properties.tmp_Nhi1
    @views tau_subcolumn[:] .= 0.0

    avx_sum_along_columns_between!(
        tau_subcolumn,
        rt.optical_properties.gas_tau[sve.gas],
        idx1,
        idx2
        )

    @turbo for i in eachindex(jac.I)

        jac.I[i] = -rt.hires_radiance.I[i] * (
            (1.0 / cosd(rt.scene.solar_zenith)) +
                (1.0 / cosd(rt.scene.observer.viewing_zenith))
        ) * tau_subcolumn[i] / scale_factor / ustrip(sve.unit, 1.0)

    end

    # Re-zero
    @views tau_subcolumn[:] .= 0.0

end

"""
$(TYPEDSIGNATURES)

Calculates the `GasLevelScalingFactorSVE` Jacobian for a `BeerLambertRTMethod` RT object
and an `UplookingGroundObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::GasLevelScalingFactorSVE,
    observer::UplookingGroundObserver
    )

    # If this gas is not part of this window, ignore
    if !haskey(rt.optical_properties.gas_tau, sve.gas)
        return
    end

    ch = create_sphericity_factors(rt.optical_properties, rt.scene)
    #sza_per_lay = create_refracted_sza(rt.optical_properties, rt.scene)

    @views jac.I[:] .= -rt.hires_radiance.I[:]

    # current value of scale factor, and apply units, then convert
    # to get value in unit of [1]
    scale_factor = sve.iterations[end] * sve.unit |> NoUnits

    idx1 = sve.start_level
    idx2 = sve.end_level - 1

    tau_subcolumn = rt.optical_properties.tmp_Nhi1
    @views tau_subcolumn[:] .= 0.0

    # total-column optical depth due to this specific gas
    #=
    avx_sum_along_columns_between!(
        tau_subcolumn,
        rt.optical_properties.gas_tau[sve.gas],
        idx1,
        idx2
        )
    =#

    @turbo for l in idx1:idx2
        #this_µ = cosd(sza_per_lay[l])
        for i in eachindex(rt.hires_radiance.I)
            #tau_subcolumn[i] += rt.optical_properties.gas_tau[sve.gas][i,l] / cosd(rt.scene.solar_zenith)
            tau_subcolumn[i] += rt.optical_properties.gas_tau[sve.gas][i,l] * ch[l]
            #tau_subcolumn[i] += rt.optical_properties.gas_tau[sve.gas][i,l] / this_µ
        end
    end


    # only single path from sun to ground observer!
    # NOTE
    # Internally, this SVE has to be a scale factor in units of [1], rather
    # than, e.g. percent, so we must back-convert to the user-specified units
    @views @. jac.I[:] = -rt.hires_radiance.I[:] * tau_subcolumn[:] / scale_factor /
        ustrip(sve.unit, 1.0) #/ cosd(rt.scene.solar_zenith)

    # Re-zero the temp array
    @views tau_subcolumn[:] .= 0.0

end


"""
$(TYPEDSIGNATURES)

Calculates the `SolarScalerPolynomialSVE` Jacobian for a `BeerLambertRTMethod` RT object.
This function should work for any kind of observer type since it scales the incident
irradiance.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::SolarScalerPolynomialSVE
    )

    # This calculation is not dependent on the observer type
    swin = rt.optical_properties.spectral_window

    if swin === sve.swin
        # Only perform the Jacobian calculation if this
        # SVE points to the RT spectral window!

        # dI / dSS = dI / dSS * dSS / dSS_coefficient =
        #          = I / SS * (ww - ww_reference)^order
        @turbo for i in eachindex(swin.ww_grid)
            jac.I[i] = rt.hires_radiance.I[i] / rt.solar_scaler[i]
        end

        if sve.coefficient_order > 0

            @turbo for i in eachindex(swin.ww_grid)
                jac.I[i] *= (swin.ww_grid[i] - swin.ww_reference) ^ sve.coefficient_order
            end

        end

    else

        jac.I[:] .= 0

    end

end

"""
$(TYPEDSIGNATURES)

Calculates the partial derivative ∂I/∂τ for a `SatelliteObserver` viewing geometry for a
`BeerLambertRTMethod` RT object.
"""
function calculate_dI_dTau(
    rt::BeerLambertRTMethod,
    observer::SatelliteObserver
)
    @debug "[RT] Calculating ∂I/∂Tau for $(rt) and $(observer)"

    return -rt.hires_radiance.I *
        (1/cosd(rt.scene.solar_zenith) + 1/cosd(observer.viewing_zenith))

end

"""
$(TYPEDSIGNATURES)

Calculates the partial derivative ∂I/∂τ for a `UplookingGroundObserver` viewing geometry
for a `BeerLambertRTMethod` RT object.
"""
function calculate_dI_dTau(
    rt::BeerLambertRTMethod,
    observer::UplookingGroundObserver
)
    @debug "[RT] Calculating ∂I/∂Tau for $(rt) and $(observer)"

    return -rt.hires_radiance.I * (1/cosd(rt.scene.solar_zenith))

end


"""
$(TYPEDSIGNATURES)

Calculates the partial derivative ∂I/∂VMR for each level for a particular **gas** of type
`GasAbsorber` for the RT object **rt** of type `BeerLambertRTMethod`. Note that this
function allocates a new vector of length `rt.scene.atmosphere.N_level`.

This function is mainly needed to produce column averaging kernels for scaling retrievals.
"""
function calculate_dI_dVMR(
    rt::BeerLambertRTMethod,
    gas::GasAbsorber
    )

    # Create output array
    dI_dVMR = zeros(
        rt.optical_properties.spectral_window.N_hires,
        rt.scene.atmosphere.N_level
        )

    # First check if the gas is present
    if !haskey(rt.optical_properties.gas_derivatives, gas)
        @warn "Gas $(gas) not present in $(rt)!"
        return dI_dVMR
    end

    # Grab ∂Tau/∂VMR for this gas, where VMR has unit 1
    dTau_dVMR = rt.optical_properties.gas_derivatives[gas]["dTau_dVMR"]

    # Grab ∂I/∂Tau
    # (this is specific to RT mode and observer!)
    dI_dTau = calculate_dI_dTau(rt)

    # Calculate ∂I/∂VMR_l where l is the level
    for l in 1:rt.scene.atmosphere.N_level
        if l == 1
            # TOA level
            # at l = 1, the index at τ refers to layers
            # ∂I/∂VMR_1 = ∂I/∂τ_1 * ∂τ_1/∂VMR1
            @views @. dI_dVMR[:,l] = dI_dTau[:] * dTau_dVMR[:,l,1]
        elseif l == rt.scene.atmosphere.N_level
            # Surace level
            @views @. dI_dVMR[:,l] = dI_dTau[:] * dTau_dVMR[:,l-1,2]
        else
            # in between
            # ∂I/∂VMRl = ∂I/∂τ_l * ∂τ_l/∂VMR_l + ∂I/∂τ_l+1 * ∂τ_l+1/∂VMR_l
            @views @. dI_dVMR[:,l] = dI_dTau[:] * (
                dTau_dVMR[:,l-1,2] + dTau_dVMR[:,l,1])
        end
    end

    return dI_dVMR

end

"""
$(TYPEDSIGNATURES)

Calculates the `GasLevelScalingFactorSVE` Jacobian for a `BeerLambertRTMethod` RT object. This
further dispatches to the correct function for the specific observer type.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::TemperatureOffsetSVE,
    )

    @views jac.I[:] .= 0

    # Dispatch to specific observer type
    calculate_rt_jacobian!(jac, rt, sve, rt.scene.observer)

end

"""
$(TYPEDSIGNATURES)

Calculates the `TemperatureOffsetSVE` Jacobian for a `BeerLambertRTMethod` RT object
and an `SatelliteObserverGroundObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::TemperatureOffsetSVE,
    observer::SatelliteObserver
)

    # Use tmp array for this calculation
    dTau_dT_sum = rt.optical_properties.tmp_Nhi1
    # Zero out
    @views dTau_dT_sum[:] .= 0.0

    # T offset affects all gases, so we iterate and sum over all
    for gas in keys(rt.optical_properties.gas_derivatives)

        # We can double-check here if a ∂τ/∂T was actually calculated
        if !haskey(rt.optical_properties.gas_derivatives[gas], "dTau_dT")
            @debug "[RT] No ∂τ/∂T found for $(gas)"
            continue
        end

        dTau_dT = rt.optical_properties.gas_derivatives[gas]["dTau_dT"]
        # Sum up contributions for layers
        avx_add_along_columns!(dTau_dT_sum, dTau_dT)


    end

    unit_fac = ustrip(rt.scene.atmosphere.temperature_unit, 1.0 * sve.unit)

    @turbo for i in eachindex(jac.I)

        jac.I[i] = -rt.hires_radiance.I[i] * (
            (1.0 / cosd(rt.scene.solar_zenith)) +
                (1.0 / cosd(rt.scene.observer.viewing_zenith))
        ) * dTau_dT_sum[i] / unit_fac

    end

end

"""
$(TYPEDSIGNATURES)

Calculates the `GasLevelScalingFactorSVE` Jacobian for a `BeerLambertRTMethod` RT object
and an `UplookingGroundObserver` observer.
"""
function calculate_rt_jacobian!(
    jac::Radiance,
    rt::BeerLambertRTMethod,
    sve::TemperatureOffsetSVE,
    observer::UplookingGroundObserver
)

    # Use tmp array for this calculation
    dTau_dT_sum = rt.optical_properties.tmp_Nhi1
    # Zero out
    @views dTau_dT_sum[:] .= 0.0

    # T offset affects all gases, so we iterate and sum over all
    for gas in keys(rt.optical_properties.gas_derivatives)

        # We can double-check here if a ∂τ/∂T was actually calculated
        if !haskey(rt.optical_properties.gas_derivatives[gas], "dTau_dT")
            @debug "[RT] No ∂τ/∂T found for $(gas)"
            continue
        end

        dTau_dT = rt.optical_properties.gas_derivatives[gas]["dTau_dT"]
        # Sum up contributions for layers
        # (and add them to dTau_dT_sum)
        avx_add_along_columns!(dTau_dT_sum, dTau_dT)


    end

    unit_fac = ustrip(rt.scene.atmosphere.temperature_unit, 1.0 * sve.unit)

    @turbo for i in eachindex(jac.I)

        jac.I[i] = -rt.hires_radiance.I[i] * (
            (1.0 / cosd(rt.scene.solar_zenith))
        ) * dTau_dT_sum[i] / unit_fac

    end

end