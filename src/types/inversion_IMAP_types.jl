"""

A type to facilitate inversions using the Iterative Maximum
A-posteriori (IMAP) solver method. See Frankenberg et al. (2005)
for details.

$(TYPEDFIELDS)

"""
struct IMAPSolver <: AbstractSolver
    "The forward model function, only takes AbstractStateVector as argument"
    forward_model::Function
    "The state vector"
    state_vector::AbstractStateVector
    "Prior covariance matrix"
    prior_covariance::AbstractMatrix
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

    function IMAPSolver(
        fm::Function,
        sv::AbstractStateVector,
        Sa::AbstractMatrix,
        max_iter::Int,
        dsigma_scale::Number,
        dispersions::Dict{<:AbstractSpectralWindow, <:AbstractDispersion},
        indices::Dict{<:AbstractSpectralWindow, <:AbstractVector},
        radiance::Radiance,
        jacobians::Dict{AbstractStateVectorElement, <:Radiance},
        measured::Dict{<:AbstractDispersion, <:AbstractVector},
        noise::Dict{<:AbstractDispersion, <:AbstractVector}
    )

        # Construct the prior covariance matrix.
        # It's convenient to have it as part of the structure,
        # since it is not changing during the inversion.

        # Note - any cross-correlations have to be added manually at this point.
        @views Sa[:,:] .= 0.0
        for i in axes(Sa, 1)
            Sa[i,i] = sv.state_vector_elements[i].prior_covariance
        end

        return new(
            fm,
            sv,
            Sa,
            max_iter,
            dsigma_scale,
            dispersions, # dispersions
            indices, # indices
            radiance, # radiance
            jacobians, # jacobians
            measured,
            noise
        )

    end
end
