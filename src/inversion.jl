include("inversion_IMAP.jl")
include("inversion_LM.jl")


"""
    Creates an alphabetically-ordered list of spectral windows
    to iterate over for consistent iteration order.


"""
function buffer_swin_iterator(s::AbstractSolver)

    swin_names = [swin.window_name for swin in keys(s.dispersions)]
    swin_order = sortperm(swin_names)

    return [swin for swin in keys(s.dispersions)][swin_order]

end


"""
Turns the Jacobians from the solver object into a freshly allocated matrix, so they
can be used in the inversion part of the algorithm.

$(TYPEDSIGNATURES)

"""
function create_K_from_solver(s::AbstractSolver)

    N_jac = length(s.jacobians)
    N_wl = map(length, values(s.indices)) |> sum

    K = zeros((N_wl, N_jac))

    for swin in buffer_swin_iterator(s)
        for (i, sve) in enumerate(s.state_vector.state_vector_elements)
            @views K[s.indices[swin], i] .= s.jacobians[sve].I[s.indices[swin]]
        end
    end

    return K

end

"""
Creates and returns the instrument noise covariance matrix

$(TYPEDSIGNATURES)
"""
function create_Se_from_solver(s::AbstractSolver; return_inverse=false)

    N_wl = map(length, values(s.indices)) |> sum
    Se_diag = zeros(N_wl)

    # Construct diagonal entries
    for swin in keys(s.indices)

        disp = s.dispersions[swin]
        Se_diag[s.indices[swin]] = (s.instrument_noise[disp][disp.index]) .^ 2

    end

    # Create diagonal matrix from those entries
    # (can be used like a regular matrix)
    Se = Diagonal(Se_diag)

    #=
    Se = Diagonal(
            vcat([
                (s.instrument_noise[s.dispersions[swin]][s.dispersions[swin].index]) .^ 2
                for swin in buffer_swin_iterator(s)
                    ]...)
                    )
    =#

    if return_inverse
        return inv(Se)
    else
        return Se
    end

end

"""
Calculates the reduced ``χ^2`` statistic for the spectral
fit within `solver`.

$(TYPEDSIGNATURES)

# Details

The reduced ``χ^2`` statistic is calculated as

```math
χ^2 = \\frac{1}{N_{\\text{spec}} - N_{\\text{sv}}}
\\sum_i^{N_{\\text{spec}}} \\frac{M_i - O_i}{\\varepsilon_i}
```

"""
function calculate_chi2(s::AbstractSolver)

    # Initialize with zero
    chi2 = Dict{AbstractSpectralWindow, Float64}()

    # We must loop over all spectral windows used for this retrieval
    for swin in keys(s.dispersions)

        chi2[swin] = 0.0

        #=
        NOTE:

        @turbo in this loop seems to speed up the computation by a factor of ~2,
        but occasionally causes the kernel to crash, probably due to a bad value
        in either the noise or the measured radiances. So for now, we omit
        @turbo in this function.
        =#

        disp = s.dispersions[swin]

        for i in eachindex(s.indices[swin])

            idx = s.indices[swin][i]
            disp_idx = disp.index[i]

            chi2[swin] += (
                (s.radiance.I[idx] - s.measured[disp][disp_idx])^2 /
                    (s.instrument_noise[disp][disp_idx])^2
            )

        end

        N_spec = length(s.indices[swin])
        N_sv  = length(s.state_vector)
        chi2[swin] /= (N_spec)

        # Some definions subtract the number of free parameters:
        # chi2[swin] /= (N_spec - N_sv - 1)

    end

    return chi2

end

"""
Returns the measured radiances of an `AbstractSolver` `s`, belonging to an
`AbstractSpectralWindow` `swin`. The optional argument `view` determines whether
to return a view onto the array (`true`), or a copy (`false`).

$(TYPEDSIGNATURES)
"""
function get_measured(
    s::AbstractSolver,
    swin::AbstractSpectralWindow;
    view=true)

    disp = s.dispersions[swin]

    if view
        return @views s.measured[disp][disp.index]
    else
        return s.measured[disp][disp.index]
    end

end


function get_measured(s::AbstractSolver)

    # Create the needed array
    N = map(length, values(s.indices)) |> sum
    measured = zeros(N);

    # Fill with the values from the measured radiances
    get_measured!(measured, s)

    return measured

end

function get_measured!(meausured, s)

    for swin in keys(s.indices)
        @views meausured[s.indices[swin]] = get_measured(s, swin)
    end

end


"""
Returns the modeled radiance currently stored in the `AbstractRTBuffer` of an
`EarthAtmosphereBuffer` `buf`. The optional argument `view` determines whether
to return a view onto the array (`view=true`), or a copy (`view=false`).

$(TYPEDSIGNATURES)
"""
function get_modeled(
    s::AbstractSolver,
    swin::AbstractSpectralWindow;
    view=true
    )

    if view
        return @views s.radiance.I[s.indices[swin]]
    else
        return s.radiance.I[s.indices[swin]]
    end
end

# Shortcut for UK spelling
get_modelled(s, swin; view=true) = get_modeled(s, swin; view)

function get_modeled(s::AbstractSolver)

    # Create the needed array
    N = map(length, values(s.indices)) |> sum
    modelled = zeros(N);

    # Fill with the values from the modelled radiances
    get_modeled!(modelled, s)

    return modelled

end

function get_modeled!(modelled, s)

    for swin in keys(s.indices)
        @views modelled[s.indices[swin]] = get_modelled(s, swin)
    end

end


function get_noise(
    s::AbstractSolver,
    swin::AbstractSpectralWindow;
    view = true
    )

    disp = s.dispersions[swin]

    if view
        return @views s.instrument_noise[disp][disp.index]
    else
        return s.instrument_noise[disp][disp.index]
    end

end

function get_noise(s::AbstractSolver)

    # Create the needed array
    N = map(length, values(s.indices)) |> sum
    noise = zeros(N);

    # Fill with the values from the instrument noise
    get_noise!(noise, s)

    return noise

end

function get_noise!(noise, s)

    for swin in keys(s.indices)
        @views noise[s.indices[swin]] = get_noise(s, swin)
    end

end


# This needs to be solved more elegantly
get_wavenumber(s) = get_wavelength(s)
get_wavenumber(s, swin) = get_wavelength(s, swin)

function get_wavelength(
    s::AbstractSolver,
    swin::AbstractSpectralWindow;
    view = true
    )

    disp = s.dispersions[swin]

    if view
        return @views disp.ww[:]
    else
        return disp.ww[:]
    end

end

function get_wavelength(s::AbstractSolver)

    # Create the needed array
    N = map(length, values(s.indices)) |> sum
    ww = zeros(N);

    for swin in keys(s.indices)
        @views ww[s.indices[swin]] = get_wavelength(s, swin)
    end

    return ww
end


"""

$(TYPEDSIGNATURES)
"""
function calculate_scale_factor_AK(
    buf::AbstractAtmosphereBuffer,
    s::AbstractSolver,
    gas::GasAbsorber,
    isrf::Dict{<:AbstractSpectralWindow, <:AbstractISRF};
    G::Union{Nothing, AbstractMatrix}=nothing
    )

    if !(any_SVE_is_type(s.state_vector, GasLevelScalingFactorSVE))
        @error "Sorry - no GasLevelScalingFactorSVE in state vector!"
        return nothing
    end
    # Grab short-hand for number of levels
    N_level = buf.scene.atmosphere.N_level


    # Which position in the state vector relates to this gas scale factor?
    SV_idx = -1
    for (idx, sve) in StateVectorIterator(s.state_vector, GasLevelScalingFactorSVE)
        if sve.gas === gas
            SV_idx = idx
        end
    end

    if SV_idx == -1
        @error "Could not find scale factor SVE for $(gas)"
        return nothing
    end

    # If no gain matrix G was supplied, calculate it here
    if isnothing(G)

        q = calculate_OE_quantities(s)
        G = q.G

    end

    #=
        Calculate ∂I/∂VMR for all spectral windows
        and apply the instrument function
    =#

    # Pre-allocate array to contain low-resolution ∂I/∂VMR
    N_wl = map(length, values(s.indices)) |> sum
    dI_dVMR_low = zeros(N_wl, N_level)


    # Calculate dI_dVMR for all spectral windows, and all levels!
    # apply the instrument function, and then store in dI_dVMR_low
    for swin in keys(s.indices)

        dI_dVMR = calculate_dI_dVMR(buf.rt[swin], gas)


        for l in 1:buf.scene.atmosphere.N_level

            x = dI_dVMR[:,l]

            # Note that here, we need the ISRF corresponding
            # to the particular band
            success = apply_isrf_to_spectrum!(
                buf.inst_buf,
                isrf[swin],
                buf.rt_buf.dispersion[swin],
                x,
                swin
            )

            # buf.rt_buf.indices[swin] makes sure we put the
            # ISRF-applied low-res spectra into the right place
            # in the array
            @views dI_dVMR_low[buf.rt_buf.indices[swin],l] =
                buf.inst_buf.low_res_output[buf.rt_buf.dispersion[swin].index]

        end

    end

    # Now calculate the pressure weighting function h
    h = create_pressure_weights(buf.scene.atmosphere)

    # Convert current VMR profile into unitless quantity
    fg_vmr = gas.vmr_levels * gas.vmr_unit .|> NoUnits

    # Scale back the VMR to its first guess value
    for (idx,sve) in StateVectorIterator(s.state_vector, GasLevelScalingFactorSVE)
        if !(sve.gas === gas)
            # Skip over non-matching gas-SVE combinations
            continue
        end

        fg_vmr[sve.start_level:sve.end_level] ./=
            get_current_value_with_unit(sve)
    end

    # Create a column vector that extracts the parts of the
    # VMR profile where the scale factor acts
    c = zeros(N_level)
    # Allocate output
    AK = zeros(N_level)

    # We have to loop over all SVEs and add up their contribution to the
    # AK. Some set-ups might include several SVEs that scale different parts
    # of the gas VMR profile.
    for (idx, sve) in StateVectorIterator(s.state_vector,
                                          GasLevelScalingFactorSVE)

        if !(sve.gas === gas)
            # Skip over non-matching gas-SVE combinations
            continue
        end

        @views c[:] .= 0
        c[sve.start_level:sve.end_level] .= 1.0

        # Unit conversion factor
        unit_fac = 1.0 * get_unit(sve) |> NoUnits

        # Calculate column scale AK as
        # G * ∂I/∂VMR * (h^T * VMR)
        # .. and a last factor to account for the unit conversion
        # between SVE unit (e.g. percent, or nothing) and unit 1.
        @debug "[INV] Adding up for $(sve) at $(idx)"
        for l in 1:N_level
            AK[l] += (G * dI_dVMR_low[:,l])[idx] * dot(h, fg_vmr .* c) * unit_fac
        end

    end

    return AK ./ h
end

"""
Print a convenient summary of the current state vector, including
the associated uncertainties.

$(TYPEDSIGNATURES)

Note that this function does not check for convergence.

"""
function print_posterior(s::AbstractSolver)

    q = calculate_OE_quantities(s)

    if !isnothing(q)
        print_posterior(q)
    else
        @error "[INV] OE quantities could not be calculated."
    end

end

# Legacy way
print_posterior(s::AbstractSolver, q::OEQuantities) = print_posterior(q)

"""
Print a convenient summary of the current state vector, including
the associated uncertainties.

$(TYPEDSIGNATURES)

Note that this function does not check for convergence.

"""
function print_posterior(q::OEQuantities)
    println("NNEEEWWW")
    sv = q.state_vector

    if !isnothing(q)

        pretty_table(
            (
                hcat(
                    get_name(sv),
                    get_unit(sv),
                    get_current_value(sv),
                    q.SV_ucert,
                    q.SV_ucert_smoothing,
                    q.SV_ucert_noise,
                    diag(q.AK)
                )
            ),
            column_labels=[
                ["Name", "Units", "Value",
                 "Uncertainty", "Uncertainty", "Uncertainty",
                 "AK"],
                ["", "", "",
                 "(total)", "(smoothing)", "(noise)",
                 ""]
            ],
            display_size=(300, 200),
            title="Posterior state vector",
        )
    else
        @error "No OEQuantities available."
    end

end

"""
Returns the number of iterations performed in an IMAPSolver
object, and also checks if all state vector elements have the
same count (as they should).

$(TYPEDSIGNATURES)

"""
function get_iteration_count(s::AbstractSolver)
    if length(s.state_vector) == 0
        @warn "This state vector has length zero (is empty)!"
        return 0
    else

        l = unique(
            map(x -> length(x.iterations),
                s.state_vector.state_vector_elements)
                   )

        if length(l) == 1

            return l[1]

        else

            @warn "Not all state vector elements have the same iteration counts!"
            return 0

        end

    end
end


function check_solver_validity(s::AbstractSolver)

    # Loop through all present spectral windows
    for swin in keys(s.indices)

        index = s.indices[swin]

        # Check if we have sensible values in the index
        if length(index) < 1
            @warn "Empty index vector in solver."
            return false
        end

        # Check for non-finites in radiance

        if check_for_not_finite(s.radiance.I[index])
            @warn "Non-finites in solver radiances."
            return false
        end

        # Check for non-finites in jacobians
        for j in eachindex(s.jacobians)
            if check_for_not_finite(s.jacobians[j].I[index])
                @warn "Non-finites in solver Jacobian number $(j)."
                return false
            end
        end

    end

    return true

end


"""
    Computes the finite-difference Jacobian for the forward model embedded in solver `s`.
    Users must supply a perturbation value `Δx` (must be unit-compatible with `sve`), such
    that (F(x + Δx) - F(x)) / Δx can be computed. `Δx` must have units.
"""
function compute_finite_difference_jacobian(
    s::AbstractSolver,
    sve::AbstractStateVectorElement,
    Δx
)

    # Check if SVE is indeed contained in the solver, otherwise something went wrong..
    is_there = false
    for solver_sve in s.state_vector.state_vector_elements
        if solver_sve === sve
            is_there = true
            break
        end
    end

    if !is_there
        error("$(sve) not found in solver state vector!")
    end

    # Create output container
    result_1 = similar(s.radiance)
    result_1[:] .= 0
    result_2 = similar(s.radiance)
    result_2[:] .= 0
    result_FD = similar(s.radiance)
    result_FD[:] .= 0

    #=
        Perturb the SVE

        IMPORTANT

        Running a forward model usually changes the atmospheric state; for example, a
        forward model could be set up to construct a new pressure grid when a surface
        pressure SVE is part of the state vector. Another example is the temperature
        offset SVE, which would add a constant offset to the temperature profile.

        Therefore, if we perturb the SVE here, and then run the forward model, the
        atmosphere will be altered, and any other

    =#

    # take a copy of the atmosphere (ignoring atmospheric elements..)

    # Add an artificial iteration to ALL state vector elements
    for _sve in s.state_vector.state_vector_elements
        push!(_sve.iterations, _sve.iterations[end])
    end

    # Compute the forward model for `x` (F(x))
    s.forward_model(s.state_vector)

    # Extract the result
    @views result_1[:,:] = s.radiance[:,:]


    @info "SVE original: $(sve)"

    # Perturb the SVE we want to analyze
    sve.iterations[end] += Δx
    @info "SVE perturbed: $(sve)"

    # Compute the forward model for `x + Δx` (F(x + Δx))
    s.forward_model(s.state_vector)

    # Extract the result and add to the numerator
    @views result_2[:,:] = s.radiance[:,:]


    # Reverse the perturbation of the SVE
    sve.iterations[end] -= Δx
    @info "SVE re-set: $(sve)"

    # Delete the artificial iteration
    for _sve in s.state_vector.state_vector_elements
        pop!(_sve.iterations)
    end
    # Divide by Δx to get [F(x + Δx) - F(x)] / Δx
    @views @. result_FD[:,:] = (result_2[:,:] - result_1[:,:]) / Δx

    # The results are F(x), F(x + Δx), [F(x + Δx) - F(x)]/Δx
    return result_1, result_2, result_FD

end