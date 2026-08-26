# [Custom instrument model](@id custom_instrument_model)

RetrievalToolbox does not provide any particular instrument model. There is, however, one function that applies instrument spectral response functions (ISRFs) to modelled radiance to obtain the instrument-level radiance (see [here](@ref ISRF_concept)). This is a very generic formalism as long as the spectral response of an instrument can be characterized via either some analytic expression or through tabulated (measured) data.

Users may want to write their own instrument model, for example if the ISRF formalism does not capture some aspect of an instrument. In this case it is important to make sure that the new custom instrument model inserts the results into the right place. Users can either study the source code of the `apply_isrf_to_spectrum!` function (`src/instrument.jl`), or read the rest of this section.

## General recipe

Let `rt` be a general RT method (see [here](@ref RT_method_types)) that contains a `hires_radiance` object which can be either a `ScalarRadiance` or a `VectorRadiance` (see [here](@ref radiance_types_intro) for details on how radiances are handled). The task of the instrument function is to take the contents of `rt`, perform the user-specific instrument model calculations, and then copy the results into the correct positions of the `RTBuffer` (see [here](@ref buffer_types)).

The correct mapping between spectral space coordinates (wavelength or wavenumber) and spectral elements of the detector is goverened by the **dispersion objects** (see [type definition](@ref dispersion_types) and a more thorough [explanation](@ref samples)), which are always linked to a specific **spectral window object**. It can be understood as a lookup which says, for example: spectral element `70` on the detector is sensitive to radiance at a central wavelength of `1.60423 µm`, and so on. Hence, any user-defined instrument model should calculate the detector-level radiance for every spectral point on the `disp.ww` vector.

!!! note "Custom instrument model function"
    For every spectral point on the instrument-level spectral grid `disp.ww`, representing the spectral samples covered by the attached spectral window, the custom instrument model must calculate the detector-level radiances. 


While the dispersion object(s) provide the relationship between spectral space and detector samples, the inversion routines require an additional look up in order to find the correct spectral samples **within the radiative transfer (RT) buffer**. This is handled by the `indices` field of the `rt` object. It is a dictionary that links an `AbstractSpectralWindow` to an array of numbers that each correspond to a positional mapping between instrument (dispersion) sample and the position of that sample in the RT buffer. Assuming all buffers, spectral windows and dispersions have been initialized correctly, and the `calculate_indices!` function has been called at some point beforehand (ideally towards the beginning of the forward model), the overall custom instrument model function should look roughly as follows:

```julia

# Assume `buf` is an `EarthAtmosphereBuffer`
rt_buf = buf.rt_buf

# Loop through every RT object in the EarthAtmosphereBuffer
for (swin, rt) in buf.rt
    

    # Grab the correct dispersion attached to the
    # spectral window `swin`
    disp = rt_buf.dispersion[swin]

    # CUSTOM INSTRUMENT MODEL CALCULATION
    # ===================================
    #
    # This is the function users MUST implement themselves
    #
    # Takes `rt.hires_radiance.I`, which is the model radiance 
    # intensity at wavelengths (or wavenumbers) `swin.ww` and
    # performs calculations to obtain instrument-level radiance
    # at `disp.ww`, the central wavelengths (or wavenumbers) at
    # the detector samples `disp.index`. Let's call the result
    # `instrument_radiance`.
    # (or use `rt.hires_radiance.S`) for the more complete 
    # Stokes vector if polarization is required.

    instrument_radiance = my_own_instrument_model(
        rt.hires_radiance.I, # high-res model radiance
        disp.ww, # detector-level wavelengths
        disp.index, # detector sample absolute indexing
        disp.detector_samples, # full detector grid
        ... # more arguments?
    )
    # Note: `instrument_radiance` must have the same length as
    # `disp.detector_samples`, or the total number of detector
    # samples for this detector.

    # Now copy over the results from `instrument_radiance` into
    # the correct postitions of the RT buffer. Use an explicit
    # loop for clarity.

    for i in 1:length(disp.index)

        # Position in the RT buffer (destination)
        rt_buf_idx = rt_buf.indices[swin][i]
        # Position in the custom calculation result (source)
        disp_idx = disp.index[i]

        # copy to buffer!
        rt_buf.radiance.I[rt_buf_idx] = instrument_radiance[disp_idx]

        # Perform unit conversion
        rt_buf.radiance.I[rt_buf_idx] *= rt.radiance_unit / rt_buf.radiance_unit

    end
            
end
```

**Note the following!** The above code snippet assumes that the user function `my_own_instrument_model` produces a new vector with the instrument-level radiance and uses absolute indexing to refer to the detector sample positions. So even if this particular spectral window only covers a small section of the detector's full range, users should ensure that their custom function produces a full vector to represent the full detector. Unused detector samples can simply be set to zero or `NaN`. Users could use their own bookkeeping here, but we **strongly recommend** to stick with the described scheme, as this would preserve compatibility with the ISRF formalism. It allows for easy switching between those two ways of implementing an instrument model.

Further, the same operation also needs to be done to convert the model level Jacobians (`rt.hires_jacobians[sve]` for state vector element `sve`) into at-instrument-Jacobians (`rt_buf.jacobians[sve]`). Note that not all Jacobians are calculated analytically when the radiances are computed, but some are handled manually. For example, the dispersion polynomial coefficient derivative must be computed via `calculate_dispersion_polynomial_jacobian!`, and the result then must be handled by the user.

## Additional details

The custom instrument model may add separate modelled effects, such as stray light. Note that users also must account for instrument Doppler shift at this stage, which is the spectral shift of observed radiance due to the relative motion between observer and radiance source (which in the case of a satellite observer, would be the ground footprint on Earth).