"""
Pretty printing for Monochromatic RT

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", rt::MonochromaticRTMethod)
    # Just print this for now
    println(io, "Monochromatic RT $(rt.model)")

end

"""
Pretty printing for Monochromatic RT

$(SIGNATURES)
"""
function show(io::IO, rt::MonochromaticRTMethod)
    # Just print this for now
    println(io, "Monochromatic RT $(rt.model)")

end


function calculate_radiances_and_jacobians!(rt::MonochromaticRTMethod)

    # Make explicit dispatch to function, depending on observer mode
    calculate_radiances_and_jacobians!(rt, rt.scene.observer)

end

function calculate_radiances_and_jacobians!(
    rt::MonochromaticRTMethod,
    observer::Union{SatelliteObserver, UplookingGroundObserver}
    )

    #=
        With monochromatic RT, we simply pass all spectral points from the
        RT object to the specific RT model interfaces
    =#

    if rt.model == :XRTM
        _calculate_radiances_and_jacobians_XRTM!(rt, observer)
    else
        throw(ArgumentError("RT model $(rt.model) is not implemented!"))
    end

end


# Default behaviour for RT jacobians
function calculate_rt_jacobian!(
    rt::MonochromaticRTMethod,
    sve::AbstractStateVectorElement
)
    @debug "No Monochromatic RT Jacobian calculated for $(sve)"
    return false
end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::AbstractStateVectorElement
)

    @warn "Not implemented for $(typeof(sve))"
    return false
end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::BRDFPolynomialSVE
    )

    #=
        Important!

        Here we decide whether to move on with the calculation!
        `sve` must point to the spectral window that is attached to this particular `rt`
        via the scene object.

    =#

    if !(sve.swin === rt.optical_properties.spectral_window)
        @views jac[:] .= 0.0
        return nothing
    end


    #=
        There is no need to dispatch to a specific observer type here.
        We will calculate ∂I/∂ρ_i as

            ∂I/∂ρ_i = ∂I/∂s * ∂s/∂ρ_i = ∂I/∂s * (ww - ww_reference)^order,

        where `s` is the surface kernel amplitude (i.e. the albedo), and `ρ_i` is the
        i'th polynomial coefficient that makes up the spectrally dependent surface albedo.

    =#

    swin = rt.optical_properties.spectral_window
    idx = rt.wfunctions_map["dI_dSurface"][1]

    # Zero out Jacobian radiance object
    @turbo jac[:] .= 0

    ∂I_∂s = rt.hires_wfunctions[idx]
    o = sve.coefficient_order

    if o == 0
        @turbo jac.S[:] .= ∂I_∂s.S[:]
    else
        @turbo for i in axes(jac, 1)
            for s in axes(jac, 2)
                jac.S[i,s] = ∂I_∂s.S[i,s] * (swin.ww_grid[i] - swin.ww_reference)^o
            end
        end
    end

end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::GasVMRProfileSVE
)

    # If this gas is not part of this window, ignore
    if !haskey(rt.optical_properties.gas_tau, sve.gas)
        @debug "Skipping RT Jacobian for $(sve.gas)"
        return
    end

    # Zero out Jacobian radiance object
    jac[:] .= 0

    # Get the unit conversion factor
    unit_fac = ustrip(sve.unit, 1.0)

    # dimensions: spectral, layer
    ω = rt.optical_properties.total_omega
    τ = rt.optical_properties.total_tau

    # dimensions: spectral, layer, (higher, lower)
    ∂τ_∂VMR = rt.optical_properties.gas_derivatives[sve.gas]["dTau_dVMR"]

    if sve.level == 1

        lay = 1
        # TOA level change
        wf_idx_tau = rt.wfunctions_map["dI_dTau"][lay]
        # This is an Radiance object
        ∂I_∂τ = rt.hires_wfunctions[wf_idx_tau]

        # Add ∂I/∂τ contributions
        @turbo for i in axes(jac, 1)
            for s in axes(jac, 2)
                jac.S[i,s] += ∂I_∂τ.S[i,s] * ∂τ_∂VMR[i,lay,1]
            end
        end
        # Only use ∂I/∂ω contributions if they exist
        if haskey(rt.wfunctions_map, "dI_dOmega")
            wf_idx_omega = rt.wfunctions_map["dI_dOmega"][lay]
            ∂I_∂ω = rt.hires_wfunctions[wf_idx_omega]
            @turbo for i in axes(jac, 1) # Spectral
                for s in axes(jac, 2) # Stokes
                    jac.S[i,s] += ∂I_∂ω.S[i,s] * (-ω[i,lay] / τ[i,lay]) * ∂τ_∂VMR[i,lay,1]
                end
            end

        end
    elseif sve.level == rt.scene.atmosphere.N_level

        # BOA level change
        lay = sve.level - 1
        wf_idx_tau = rt.wfunctions_map["dI_dTau"][lay]
        ∂I_∂τ = rt.hires_wfunctions[wf_idx_tau]

        # Add ∂I/∂τ contributions
        @turbo for i in axes(jac, 1) # Spectral
            for s in axes(jac, 2) # Stokes
                jac.S[i,s] += ∂I_∂τ.S[i,s] * ∂τ_∂VMR[i,lay,2]
            end
        end

        # Only use ∂I/∂ω contributions if they exist
        if haskey(rt.wfunctions_map, "dI_dOmega")
            wf_idx_omega = rt.wfunctions_map["dI_dOmega"][lay]
            ∂I_∂ω = rt.hires_wfunctions[wf_idx_omega]

            @turbo for i in axes(jac, 1) # Spectral
                for s in axes(jac, 2) # Stokes
                    jac.S[i,s] += ∂I_∂ω.S[i,s] * (-ω[i,lay] / τ[i,lay]) * ∂τ_∂VMR[i,lay,2]
                end
            end
        end
    else
        # Everything in between
        #=
            Layers start at index 1, levels at index 1 as well,
            so for a given level (i) where i > 1, the layer above is (i-1) and
            the layer below is (i).

            Further - the change in VMR at level (l), incorporates the ∂τ/∂VMR
            for two layers. The change of τ in the layer above (l-1) w.r.t. the
            VMR change in the level below said layer; and the change of τ in layer
            below (l) w.r.t. the VMR change in the level above.

        =#
        wf_idx_tau_up = rt.wfunctions_map["dI_dTau"][sve.level - 1]
        wf_idx_tau_lo = rt.wfunctions_map["dI_dTau"][sve.level]
        l_up = sve.level - 1
        l_lo = sve.level
        ∂I_∂τ_up = rt.hires_wfunctions[wf_idx_tau_up]
        ∂I_∂τ_lo = rt.hires_wfunctions[wf_idx_tau_lo]

        # ∂I/∂τ * ∂τ/∂VMR, and account for level-layer chain rule
        @turbo for i in axes(jac, 1) # Spectral
            for s in axes(jac, 2) # Stokes
                jac.S[i,s] += ∂I_∂τ_up.S[i,s] * ∂τ_∂VMR[i,l_up,1]
                jac.S[i,s] += ∂I_∂τ_lo.S[i,s] * ∂τ_∂VMR[i,l_lo,2]
            end
        end

        # Only use ∂I/∂ω contributions if they exist
        # (+ ∂I/∂ω * ∂ω/∂VMR = + ∂I/∂ω * (-ω/τ) * ∂τ/∂VMR)
        if haskey(rt.wfunctions_map, "dI_dOmega")

            wf_idx_omega_up = rt.wfunctions_map["dI_dOmega"][sve.level - 1]
            wf_idx_omega_lo = rt.wfunctions_map["dI_dOmega"][sve.level]
            ∂I_∂ω_up = rt.hires_wfunctions[wf_idx_omega_up]
            ∂I_∂ω_lo = rt.hires_wfunctions[wf_idx_omega_lo]

            @turbo for i in axes(jac, 1) # Spectral
                for s in axes(jac, 2) # Stokes
                    jac.S[i,s] += ∂I_∂ω_up.S[i,s] * (- ω[i,l_up] / τ[i,l_up]) *
                        ∂τ_∂VMR[i,l_up,1]
                    jac.S[i,s] += ∂I_∂ω_lo.S[i,s] * (- ω[i,l_lo] / τ[i,l_lo]) *
                        ∂τ_∂VMR[i,l_lo,2]
                end
            end
        end
    end

    # Make sure to apply the unit factor to all
    @turbo jac[:] ./= unit_fac

end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::GasLevelScalingFactorSVE
    )

    # If this gas is not part of this window, ignore
    if !haskey(rt.optical_properties.gas_tau, sve.gas)
        @debug "Skipping RT Jacobian for $(sve.gas)"
        return
    end

    # current value of scale factor, and apply units
    scale_factor = sve.iterations[end] * sve.unit |> NoUnits
    # Factor due to unit
    unit_factor = ustrip(sve.unit, 1.0)

    # Zero out, since we add up layer contributions
    jac[:] .= 0

    # Layer summation start and end
    idx1 = sve.start_level
    idx2 = sve.end_level - 1
    τ_gas = rt.optical_properties.gas_tau[sve.gas]
    ω = rt.optical_properties.total_omega
    τ = rt.optical_properties.total_tau

    for l in idx1:idx2

        # Grab the index to access the appropriate RT weighting function
        wf_idx_tau = rt.wfunctions_map["dI_dTau"][l]
        # This is a Radiance
        ∂I_∂τ = rt.hires_wfunctions[wf_idx_tau]
        if haskey(rt.wfunctions_map, "dI_dOmega")
            wf_idx_omega = rt.wfunctions_map["dI_dOmega"][l]
            ∂I_∂ω = rt.hires_wfunctions[wf_idx_omega]
        end

        # Remember: VMR ~ c * τ_orig (c is the scale factor)
        # so that ∂VMR/∂c = τ/c, and we must calculate
        # ∂I/∂VMR = ∂I/∂τ * ∂τ/∂VMR + ∂I/∂ω * ∂ω/∂VMR

        @turbo for i in axes(jac, 1) # Spectral
            for s in axes(jac, 2) # Stokes
                jac.S[i,s] += ∂I_∂τ.S[i,s] * τ_gas[i,l] / scale_factor / unit_factor
            end
        end

        if haskey(rt.wfunctions_map, "dI_dOmega")
            @turbo for i in axes(jac, 1) # Spectral
                for s in axes(jac, 2) # Stokes
                    jac.S[i,s] += ∂I_∂ω.S[i,s] * (-ω[i,l] / τ[i,l]) *
                        τ_gas[i,l] / scale_factor / unit_factor
                end
            end
        end

    end

end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::SurfacePressureSVE
    )

    # Need dI_dTau at least ..
    if !("dI_dTau" in keys(rt.wfunctions_map))
        @error "dI_dTau weighting functions not found!"
    end


    # What gases do we have with dTau/dpsurf?
    gases = keys(rt.optical_properties.gas_derivatives)

    unit_factor = ustrip(rt.scene.atmosphere.pressure_unit, 1.0 * get_unit(sve))

    wf_idx_tau = rt.wfunctions_map["dI_dTau"][rt.scene.atmosphere.N_layer]
    τ = @view rt.optical_properties.total_tau[:,rt.scene.atmosphere.N_layer]
    ∂I_∂τ = rt.hires_wfunctions[wf_idx_tau]

    # We might or might not need partial derivative w.r.t. ω
    if "dI_dOmega" in keys(rt.wfunctions_map)
        wf_idx_omega = rt.wfunctions_map["dI_dOmega"][rt.scene.atmosphere.N_layer]
        ∂I_∂ω = rt.hires_wfunctions[wf_idx_omega]
        ω = @view rt.optical_properties.total_omega[:,rt.scene.atmosphere.N_layer]
    else
        ∂I_∂ω = nothing
        ω = nothing
    end

    # Zero-out; we don't need any previous results
    jac[:] .= 0

    #=
        (the sum runs over all gases)
        ∂I/∂psurf = ∑ ∂I/∂τ * ∂τ/∂psurf + ∂I/∂ω * ∂ω/∂psurf

        where

        ∂ω/∂psurf = ∑ -ω/τ * ∂τ/∂psurf
    =#

    for gas in gases

        ∂τ_∂psurf = rt.optical_properties.gas_derivatives[gas]["dTau_dpsurf"]


        # Dispatch to specialized function that calculates ∂I/∂psurf
        # (needed for performance as objects in this function scope are not type stable)

        # Note that this is adding contributions to `jac`

        _calculate_rt_jacobian_psurf!(
            jac,
            ∂I_∂τ,
            ∂τ_∂psurf,
            ∂I_∂ω,
            ω,
            τ,
            unit_factor
        )

    end
end

function _calculate_rt_jacobian_psurf!(
    jac::Radiance,
    ∂I_∂τ::Radiance,
    ∂τ_∂psurf::Vector,
    ∂I_∂ω, # Partial deriv. w.r.t. ω [wavelength, stokes]
    ω, # Single scatter albedo for surface layer [wavelength]
    τ, # Total optical depth for surface layer [wavelength]
    unit_factor
)

    @turbo for i in axes(jac, 1)
        for s in axes(jac, 2)
            jac.S[i,s] += ∂I_∂τ.S[i,s] * ∂τ_∂psurf[i] * unit_factor
        end
    end
    # ∂I/∂ω * ∂ω/∂psurf (if needed, where ∂ω/∂psurf = -ω/τ * ∂τ/∂psurf)
    if !isnothing(ω)
        @turbo for i in axes(jac, 1)
            for s in axes(jac, 2)
                jac.S[i,s] += ∂I_∂ω.S[i,s] * (-ω[i] / τ[i]) * ∂τ_∂psurf[i] * unit_factor
            end
        end
    end

end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::TemperatureOffsetSVE
    )

    # Zero out, since we add up layer contributions
    jac[:] .= 0

    # Calculate unit factor to translate between atmosphere T units and
    # statevector T units
    unit_fac = ustrip(rt.scene.atmosphere.temperature_unit, 1.0 * sve.unit)
    # Short-cuts
    τ = rt.optical_properties.total_tau
    ω = rt.optical_properties.total_omega

    # Must do this calculation for every gas present
    for gas in keys(rt.optical_properties.gas_derivatives)

        # We can double-check here if a ∂τ/∂T was actually calculated for this gas
        if !haskey(rt.optical_properties.gas_derivatives[gas], "dTau_dT")
            @warn "No ∂τ/∂T found for $(gas)"
            continue
        end

        for l in 1:rt.scene.atmosphere.N_layer

            wf_idx_tau = rt.wfunctions_map["dI_dTau"][l]
            # [spectral, stokes]
            ∂I_∂τ = rt.hires_wfunctions[wf_idx_tau]
            # [spectral, layer]
            ∂τ_∂T = rt.optical_properties.gas_derivatives[gas]["dTau_dT"]

            @turbo for i in axes(jac, 1)
                for s in axes(jac, 2)
                    jac.S[i,s] += ∂I_∂τ.S[i,s] * ∂τ_∂T[i,l] * unit_fac
                end
            end

            if haskey(rt.wfunctions_map, "dI_dOmega")
                wf_idx_omega = rt.wfunctions_map["dI_dOmega"][l]
                ∂I_∂ω = rt.hires_wfunctions[wf_idx_omega]

                @turbo for i in axes(jac, 1) # Spectral
                    for s in axes(jac, 2) # Stokes
                        jac.S[i,s] += ∂I_∂ω.S[i,s] * (-ω[i,l] / τ[i,l]) *
                            ∂τ_∂T[i,l] * unit_fac
                    end
                end
            end

        end
    end


end


function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::AerosolOpticalDepthSVE
)


    if !haskey(rt.wfunctions_map, sve)
        error("RT weighting function mapper is not aware of $(sve)")
    end

    wf_idx = rt.wfunctions_map[sve][1]
    this_wf = rt.hires_wfunctions[wf_idx]
    # Zero out
    @views jac[:] .= 0

    # Simply copy over, both `jac` and `this_wf` are the same type of
    # either ScalarRadiance or VectorRadiance, so we can simply do a
    # copy operation that broadcasts over all Stokes components.
    @turbo jac[:] .= this_wf[:]
    # This is now ∂I/∂AOD

    #=
        If we retrieve log(AOD) rather than AOD, we must correct:

        ∂I/∂log(AOD) = ∂I/∂AOD * AOD

        Of course if the state vector element uses log(AOD), then the value in
        sve.iterations is log(AOD), and we must exponentiate to get the AOD itself.
    =#
    if sve.log

        this_aod = get_current_value_with_unit(sve) |> NoUnits
        @turbo jac[:,:] .*= exp(this_aod)

    end

end

function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::AerosolHeightSVE
)


    if !haskey(rt.wfunctions_map, sve)
        error("RT weighting function mapper is not aware of $(sve)")
    end

    wf_idx = rt.wfunctions_map[sve][1]
    this_wf = rt.hires_wfunctions[wf_idx]
    # Zero out
    @views jac[:] .= 0

    # Simply copy over, both `jac` and `this_wf` are the same type of
    # either ScalarRadiance or VectorRadiance, so we can simply do a
    # copy operation that broadcasts over all Stokes components.
    @turbo jac[:] .= this_wf[:]
    # This is now ∂I/∂Height

    #=
        If we retrieve log(Height) rather than Height, we must correct:

        ∂I/∂log(Height) = ∂I/Height * Height

        Of course if the state vector element uses log(Height), then the value in
        sve.iterations is log(Height), and we must exponentiate to get the Height itself.
    =#
    if sve.log

        this_height = get_current_value_with_unit(sve) |> NoUnits
        @turbo jac[:] .*= exp(this_height)

    end

end


function calculate_rt_jacobian!(
    jac::Radiance,
    rt::MonochromaticRTMethod,
    sve::AerosolWidthSVE
)


    if !haskey(rt.wfunctions_map, sve)
        error("RT weighting function mapper is not aware of $(sve)")
    end

    wf_idx = rt.wfunctions_map[sve][1]
    this_wf = rt.hires_wfunctions[wf_idx]
    # Zero out
    @views jac[:] .= 0

    # Simply copy over, both `jac` and `this_wf` are the same type of
    # either ScalarRadiance or VectorRadiance, so we can simply do a
    # copy operation that broadcasts over all Stokes components.
    @turbo jac[:] .= this_wf[:]
    # This is now ∂I/∂Width

    #=
        If we retrieve log(Width) rather than Width, we must correct:

        ∂I/∂log(Width) = ∂I/Width * Width

        Of course if the state vector element uses log(Width), then the value in
        sve.iterations is log(Width), and we must exponentiate to get the Width itself.
    =#
    if sve.log

        this_width = get_current_value_with_unit(sve) |> NoUnits
        @turbo jac[:] .*= exp(this_width)

    end

end