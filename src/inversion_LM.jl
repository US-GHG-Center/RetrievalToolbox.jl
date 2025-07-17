"""
For an LMSolver type, calculate the next iteration, compute the
relevant quantities and update the state vector.

$(SIGNATURES)
"""
function next_iteration!(
    s::LMSolver;
    fm_kwargs=()
    )

    fm_success = s.forward_model(s.state_vector; fm_kwargs...)

    if !fm_success
        @debug "[INV] Forward model not successfully run."
        return false
    end

    if !check_solver_validity(s)
        @debug "[INV] Invalid results found in solver object."
        return false
    end

    # Create Jacobian matrix
    K = create_K_from_solver(s)

    if any(.!isfinite.(K))
        @info "[INV] Non-finites in Jacobian matrix!"
        return false
    end

    # Inverse of error covariance
    Se_inv = create_Se_from_solver(s, return_inverse=true)

    # Compile the measurement vector and the vector of radiances!
    modeled = get_modeled(s)
    measured = get_measured(s)

    # Invert regular and scaled prior covariance matrices
    Sa_inv = inv(s.prior_covariance)

    # Current and last state vector as "proper" (unit-less) vector
    sv_ap = get_prior_value(s.state_vector)
    sv_this = get_current_value(s.state_vector)

    # Change in state vector. Solve this with the ACOS-type SVD approach

    Δx = solve_LM_scaled(
        s.gamma,
        s.prior_covariance,
        K,
        Se_inv,
        modeled,
        measured,
        sv_ap,
        sv_this,
        rcond=1e-12
    )

    # Calculate cost function forecast using Jacobian
    modeled_fc = modeled + K * Δx
    c_fc = (measured - modeled_fc)' * Se_inv * (measured - modeled_fc) +
    (sv_ap - sv_this)' * Sa_inv * (sv_ap - sv_this)
    push!(s.costfunction_forecast, c_fc)

    # Calculate cost function using actual values
    c = (measured - modeled)' * Se_inv * (measured - modeled) +
        (sv_ap - sv_this)' * Sa_inv * (sv_ap - sv_this)
    push!(s.costfunction_history, c)

    # Assess whether reduction of cost function matches with the linear forecast. This
    # obviously only works from the second iteration onwards.
    # See ACOS L2 ATBD Section 3.5 (Inverse Method)
    if get_iteration_count(s) > 1

        # Note that we want the costfunction forecast for this iteration, as computed
        # during the last iteration, hence why we grab the penultimate value.
        # (the last value `s.costfunction_forecast[end]` contains the forecast
        #  produced at this iteration, valid for the next one..)
        R = (s.costfunction_history[end - 1] - s.costfunction_history[end]) /
            (s.costfunction_history[end - 1] - s.costfunction_forecast[end - 1])

        # We only do divergent steps for γ that are not too large.
        if (R < 0.0001) & (s.gamma > 1e-4) & (s.gamma < 1e4)

            # Non-linear regime, and divergent step -> iincrease γ
            s.gamma *= 10
            s.divergent_steps += 1
            @info "Divergent step! R: $(R), γ: $(s.gamma)"

            # Remove the constfunction values that we pushed in, as well as the
            # values we stored during the evaluation of the last iteration, which
            # we will evaluate again.
            pop!(s.costfunction_forecast)
            pop!(s.costfunction_forecast)
            pop!(s.costfunction_history)

            for sve in s.state_vector.state_vector_elements
                # Remove the last iteration from all state vector elements!
                # ("last" meaning the one we just evaluated)
                pop!(sve.iterations)
            end

            # Create new SV update, using new γ, but radiances and Jacobians
            # from the last valid iteration. This way, we don't need to re-evaluate
            # the forward model at the last valid step.

            # `sv_this` is now the SV at the last valid step
            # (since we removed this current iteration)
            sv_this = get_current_value(s.state_vector)
            # Grab the Jacobian matrix from the last valid step
            K_last = s.last_valid_jacobian[1:size(K,1), 1:size(K,2)]

            N = length(modeled)
            modeled_last = zeros(N)
            @views modeled_last[1:N] = s.last_valid_radiance.I[1:N]

            # This is the new SV update: K and radiances from last valid iteration,
            # but new γ!

            Δx = solve_LM_scaled(
                s.gamma,
                s.prior_covariance,
                K_last,
                Se_inv,
                modeled_last,
                measured,
                sv_ap,
                sv_this,
                rcond=1e-12
            )

            if any(isnan.(Δx))
                @error "NaNs in state vector update."
                return false
            end

            # Add new values to the prior state vector obtain (i+1)
            for (i, sve) in enumerate(s.state_vector.state_vector_elements)
                # Push back into state vector
                push!(sve.iterations, sve.iterations[end] + Δx[i])
            end

            # We must also update the costfunction forecast using this new state vector
            # update.
            modeled_fc = modeled_last + K_last * Δx
            c_fc = (measured - modeled_fc)' * Se_inv * (measured - modeled_fc) +
                (sv_ap - sv_this)' * Sa_inv * (sv_ap - sv_this)
            push!(s.costfunction_forecast, c_fc)

            return true

        elseif (R > 0.0001) & (R < 0.25) & (s.gamma < 1e4)
            # Non-linear regime -> Increase γ to make step size smaller
            s.gamma *= 10
            @info "R: $(R), γ: $(s.gamma)"
        elseif (R > 0.25) & (R < 0.75)
            # No change
            @info "R: $(R), γ: $(s.gamma)"
        elseif R > 0.75
            # Linear regime -> Reduce γ to move closer to Gauss-Newton
            s.gamma /= 2
            @info "R: $(R), γ: $(s.gamma)"
        else
            @info "No update to γ .. just keep iterating .."
            @info "R: $(R), γ: $(s.gamma)"
        end

    end # End if (iteration > 1)


    # This is a valid radiance now, so we can store the radiance into
    # `last_valid_radiance`
    @views s.last_valid_radiance[:] .= s.radiance[:]
    # Same for the Jacobians matrix..
    @views s.last_valid_jacobian[1:size(K,1), 1:size(K,2)] .= K

    # Add new values to the prior state vector obtain (i+1)
    for (i, sve) in enumerate(s.state_vector.state_vector_elements)
        # Push back into state vector
        push!(sve.iterations, sve.iterations[end] + Δx[i])
    end

    if any(isnan.(Δx))
        @error "NaNs in state vector update."
        return false
    end

    return true

end

"""
"""
function solve_LM_scaled(
    γ::Real,
    Sa::AbstractMatrix,
    K::AbstractMatrix,
    Se_inv::AbstractMatrix,
    Fx::Vector,
    y::Vector,
    xa::Vector,
    xi::Vector;
    rcond=1e-12
    )


    sigma_ap = sqrt.(diag(Sa))
    N = Diagonal(sigma_ap);
    Sa_inv = inv(Sa)

    # Compute the left-hand side
    lhs = N' * ((1 + γ) * Sa_inv + K' * Se_inv * K) * N

    # Compute the right-hand side
    rhs = N' * (K' * Se_inv * (y .- Fx) + Sa_inv * (xa .- xi))

    # Perform SVD on the LHS matrix
    U, S, V = svd(lhs, full=true)

    # The solution to Ax = b (or lhs * x = rhs) is then V' * diag(S^-1) * (U' * b)
    # But with SVD, we can trim very low singular values.
    S_inv = 1 ./ S
    for i in 2:length(S_inv)
        if S[i] <= S[1] * rcond
            S_inv[i] = 0
        end
    end

    # Calculate the state vector increment in scaled-spaces
    Δx_scaled = V * Diagonal(S_inv) * (U' * rhs)
    # Re-scale to get back to normal and return
    return N * Δx_scaled

end


"""
Checks for convergence of an IMAPSolver type.

$(TYPEDSIGNATURES)
"""
function check_convergence(s::LMSolver; verbose=false)

    # Need at least one iteration to check for convergence
    if length(s.state_vector.state_vector_elements[1].iterations) == 1
        return false
    end

    if !check_solver_validity(s)
        return false
    end

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
        @info "Δσ² = $(dsigma_square) ($(length(s.state_vector) * s.dsigma_scale))"
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
function calculate_OE_quantities(s::LMSolver)

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



