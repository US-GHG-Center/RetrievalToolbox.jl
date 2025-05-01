"""
$(TYPEDSIGNATURES)

Updates an EarthScene object's surfaces according to the state vector (current) values.
"""
function surfaces_statevector_update!(
    scene::EarthScene,
    SV::AbstractStateVector
    )

    found_any = false

    # Update Lambertian surfaces, if applicable
    for (idx, sve) in StateVectorIterator(SV, SurfaceAlbedoPolynomialSVE)

        if sve.swin in keys(scene.surfaces)
            # Only perform an update if the surface is compatible
            if scene.surfaces[sve.swin] isa LambertianPolynomialSurface

                o = sve.coefficient_order + 1
                scene.surfaces[sve.swin].coefficients[o] = get_current_value(sve)
                found_any = true

            end
        end
    end

    # Update BRDF kernels, if applicable
    for (idx, sve) in StateVectorIterator(SV, BRDFPolynomialSVE)
        # Is this the correct spectral window and is the surface a BRDF kernel?
        if (sve.swin in keys(scene.surfaces)) & 
            (scene.surfaces[sve.swin] isa BRDFSurface)

            # Loop through the kernels
            for kernel in get_surface(scene, sve.swin).kernels
                # Check if it has the right type
                if kernel isa sve.BRDF_type
                   # Update values!
                   o = sve.coefficient_order + 1
                   kernel.coefficients[o] = get_current_value(sve)
                   found_any = true
                end
            end

        end
    end

    if !found_any
        @warn "Called surfaces_statevector_update! but no " *
            "suitable updates were performed!"
    end

end


"""
$(TYPEDSIGNATURES)

Returns the surface of this `scene` that is attached to the spectral window `swin`.
"""
function get_surface(
    scene::EarthScene,
    swin::AbstractSpectralWindow
    )

    if swin isa SpectralWindow
        return scene.surfaces[swin]
    elseif swin isa BinnedSpectralWindow
        return scene.surfaces[swin.original_window]
    end

end