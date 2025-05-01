"""

A type to facilitate inversions using the Levenberg-Marquadt solver method as used in
the ACOS algorithm.

$(TYPEDFIELDS)

"""
mutable struct LMSolver <: AbstractSolver
    "The forward model function, only takes AbstractStateVector as argument"
    forward_model::Function
    "The state vector"
    state_vector::AbstractStateVector
    "Prior covariance matrix"
    prior_covariance::AbstractMatrix
    "Levenberg-Marquadt parameter γ"
    gamma::AbstractFloat
    "Number of allowed maximal iterations"
    max_iterations::Int
    "Delta-sigma scale value used for convergence checking"
    dsigma_scale::AbstractFloat
    "Dispersion objects needed to link forward model output to measurement spectral samples"
    dispersions::Dict{<:AbstractSpectralWindow, <:AbstractDispersion}
    "Indices refer to where in the buffer we store radiances/jacobians"
    indices::Dict{<:AbstractSpectralWindow, <:AbstractVector}
    "Radiance values (current iteration only)"
    radiance::Radiance
    "Jacobian values (current iteration only)"
    jacobians::Dict{<:AbstractStateVectorElement, <:Radiance}
    "Measured radiance (full-detector range)"
    measured::Dict{<:AbstractDispersion, <:AbstractVector}
    "Measurement noise belonging to the measured radiance (full-detector range)"
    instrument_noise::Dict{<:AbstractDispersion, <:AbstractVector}
    "Cost function history at iteration [i]"
    costfunction_history::Vector{<:AbstractFloat}
    "Cost function forecast for iteration [i+1]"
    costfunction_forecast::Vector{<:AbstractFloat}
    "Last valid forward model evaluation"
    last_valid_radiance::Radiance
    "Last valid Jacobian matrix"
    last_valid_jacobian::Matrix
    "Number of divergent steps so far"
    divergent_steps::Integer

    function LMSolver(
        fm::Function,
        sv::AbstractStateVector,
        Sa::AbstractMatrix,
        gamma::AbstractFloat,
        max_iter::Int,
        dsigma_scale::Number,
        dispersions::Dict{<:AbstractSpectralWindow, <:AbstractDispersion},
        indices::Dict{<:AbstractSpectralWindow, <:AbstractVector},
        radiance::Radiance,
        jacobians::Dict{<:AbstractStateVectorElement, <:Radiance},
        measured::Dict{<:AbstractDispersion, <:AbstractVector},
        noise::Dict{<:AbstractDispersion, <:AbstractVector}
    )

        # Construct the prior covariance matrix.
        # It's convenient to have it as part of the structure,
        # since it is not changing during the inversion.

        # Note - any cross-correlations have to be added manually later on.
        @views Sa[:,:] .= 0.0
        for i in axes(Sa, 1)
            Sa[i,i] = sv.state_vector_elements[i].prior_covariance
        end

        # Create a copy of the model radiance vector for use to be the
        # last valid radiance (needed for divergent, rejected steps).
        RadType = typeof(radiance).name.wrapper
        last_valid_radiance = RadType(eltype(radiance), length(radiance))
        @views last_valid_radiance[:] .= 0

        # Create an empty matrix to be used to keep the last valid Jacobian
        Nsv = length(sv)
        last_valid_jacobian = zeros(eltype(radiance), length(radiance), Nsv)

        return new(
            fm,
            sv,
            Sa,
            gamma,
            max_iter,
            dsigma_scale,
            dispersions, # dispersions
            indices, # indices
            radiance, # radiance
            jacobians, # jacobians
            measured,
            noise,
            Float64[],
            Float64[],
            last_valid_radiance,
            last_valid_jacobian,
            0
        )

    end
end
