"""
For a IMAPSolver type, calculate the next iteration, compute the
relevant quantities and update the state vector.

$(SIGNATURES)
"""
function next_iteration!(s::IMAPSolver; fm_kwargs=())

    # Evaluate forward model with current SV
    fm_success = s.forward_model(s.state_vector; fm_kwargs...)

    if !fm_success
        @debug "[INV] Forward model not successfully run."
        return false
    else
        @debug "[INV] Forward model successfully run!"
    end

    if !check_solver_validity(s)
        @debug "[INV] Invalid results found in solver object."
        return false
    end

    # Create matrices
    # TODO:
    # write some nice inline-documentation that renders
    # nicely onto a webpage?


    # This is a common issue, Sa is not invertible, so let's make checks

    Sa_inv = inv(s.prior_covariance)
    Se_inv = create_Se_from_solver(s, return_inverse=true)

    K = create_K_from_solver(s)
    Shat_inv = Sa_inv + K' * (Se_inv * K)
    Shat = inv(Shat_inv)

    # Current and last state vector as "proper" vector
    sv_ap = get_prior_value(s.state_vector)
    #sv_first = get_first_guess(s.state_vector)
    sv_this = get_current_value(s.state_vector)

    # Compute change in state vector
    # (Rodgers 2000, eq. 5.9 on page 85)

    # Compile the measurement vector for the windows
    # Note that this requires the dispersion!
    N = map(length, values(s.indices)) |> sum

    #= Important note!
        At the moment, we assume that all of our retrievals only
        care about intensity in the end, meaning that our detectors count
        photons, and not measure the polarization state. As such, the insrument's
        polarization sensitivity must be accounted for by the user, and collapsed
        into a ScalarRadiance object in the end.
    =#
    measured = zeros(N);
    modelled = zeros(N);

    for (swin, indices) in s.indices
        @views measured[indices] = get_measured(s, swin)
        @views modelled[s.indices[swin]] = get_modeled(s, swin)
    end

    delta_sv = Shat * (K' * Se_inv) *
        (measured - modelled + (K * (sv_this - sv_ap)))

    # Add new values to the prior state vector obtain (i+1)
    for (i, sve) in enumerate(s.state_vector.state_vector_elements)

        # Push back into state vector
        push!(sve.iterations, sve.prior_value + delta_sv[i])

    end

    if any(isnan.(delta_sv))
        @error "[INV] NaNs in state vector update."
        return false
    end

    return true

end


"""
Checks for convergence of an IMAPSolver type.

$(TYPEDSIGNATURES)
"""
function check_convergence(s::IMAPSolver; verbose=false)

    # Need at least one iteration to check for convergence
    if length(s.state_vector.state_vector_elements[1].iterations) == 1
        return false
    end

    if !check_solver_validity(s)
        return false
    end

    # L1B index needed to match detector radiance
    idx = s.indices
    Sa_inv = inv(s.prior_covariance)
    Se_inv = create_Se_from_solver(s, return_inverse=true)

    # Create jacobian matrix
    K = create_K_from_solver(s)
    # Calculate inverse of posterior covariance matrix
    Shat_inv = Sa_inv + K' * (Se_inv * K)

    # Change in state vector compared to last step
    delta_sv = (
        [x.iterations[end] - x.iterations[end-1] for
         x in s.state_vector.state_vector_elements]
    )

    # dsigma_square: change in state vector as a fraction of posterior
    # covariance.
    dsigma_square = dot(delta_sv, Shat_inv * delta_sv)

    if verbose
        @info "[INV] Δσ² = $(dsigma_square) ($(length(s.state_vector) * s.dsigma_scale))"
    end
     # Scale by user-defined value (dsigma_scale) and see if it meets
    # our convergence criterion.
    return dsigma_square < length(s.state_vector) * s.dsigma_scale

end

"""
Given an IMAPSolver, this function calculates the usual
post-retrieval error analytics, such as the AK matrix. All of this
is taken from Rodgers (2000) Chapter 3.2, and Frankenberg et al. (2005).
https://doi.org/10.5194/acp-5-9-2005

$(SIGNATURES)

There are a decent number of array allocations inside this function,
but given the (generally) smaller size, and the fact that they only
need to be created once the retrieval has converged, the performance
and memory penalties are small.

"""
function calculate_OE_quantities(s::IMAPSolver)

    if length(s.indices) == 0
        return nothing
    end

    if !check_solver_validity(s)
        return nothing
    end

    # Grab inverse of prior covariance needed for Kᵀ Se⁻¹ K + Sa⁻¹
    Sa_inv = inv(s.prior_covariance)
    # Instrument noise covariance
    Se = create_Se_from_solver(s)
    Se_inv = inv(Se)
    # Create jacobian matrix
    K = create_K_from_solver(s)
    # Calculate inverse of posterior covariance matrix
    Shat_inv = Sa_inv + K' * (Se_inv * K)
    Shat = inv(Shat_inv)
    # Calculate gain matrix
    G = Shat * K' * Se_inv
    # Calculate averaging kernel matrix
    AK = G * K

    # State vector "error" (uncertainty) taken straight
    # from the diagonal elements of Shat.
    SV_ucert = sqrt.(diag(Shat))

    # Uncertainty due to (random) instrument noise
    SV_ucert_noise = sqrt.(diag(G * Se * G'))

    # Uncertainty due to smoothing
    eye = I(length(s.state_vector)) # Unit matrix, size of the state vector
    SV_ucert_smoothing = sqrt.(abs.(diag((AK - eye) *
        s.prior_covariance * (AK - eye)')))

    # Alternative calculation
    # (maybe less prone to linear algebra problems with near-zero entries?)
    # SV_ucert_smoothing = @. sqrt(SV_ucert^2 - SV_ucert_noise^2)

    # They add up in quadrature, such that
    # @. SV_ucert == sqrt(SV_ucert_noise^2 + SV_ucert_smoothing^2)

    # Build type and return
    return OEQuantities(
        s.state_vector,
        K,
        Se,
        s.prior_covariance,
        Shat,
        G,
        AK,
        SV_ucert,
        SV_ucert_noise,
        SV_ucert_smoothing
    )

end



