"""
A type to facilitate inversions using the Iterative Maximum A-posteriori (IMAP) solver
method. See [Frankenberg et al. (2005)](https://doi.org/10.5194/acp-5-9-2005) for details.
This function includes a constructor which populates the supplied prior covariance matrix
with the prior covariance values of the state vector elements. Any off-diagonal covarainces
must be added **manually**, and **after** the solver object has been instantiated.

* `forward_model::Function`: The forward model function, only takes AbstractStateVector as argument
* `state_vector::AbstractStateVector`: The state vector
* `prior_covariance::AbstractMatrix`: Prior covariance matrix
* `max_iterations::Integer`: Number of allowed maximal iterations
* `dsigma_scale::AbstractFloat`: Delta-sigma scale value used for convergence checking
* `dispersions::Dict{<:AbstractSpectralWindow, <:AbstractDispersion}`: Dispersion objects needed to link forward model output to measurement spectral samples
* `indices::Dict{<:AbstractSpectralWindow, <:AbstractVector}`: Indices refer to where in the buffer we store radiances/jacobians
* `radiance::Radiance`: Radiance values (current iteration only)
* `jacobians::Dict{<:AbstractStateVectorElement, <:Radiance}`: Jacobian values (current iteration only)
* `measured::Dict{<:AbstractDispersion, <:AbstractVector}`: Measured radiance (full-detector range)
* `instrument_noise::Dict{<:AbstractDispersion, <:AbstractVector}`: Measurement noise belonging to the measured radiance (full-detector range)

# Forward model function

In any `IMAPSolver` object, a forward model function has be supplied that returns a `Bool`
to signify whether the forward model execution was successful. By definition, the first
and only non-keyword argument must be an `AbstractStateVector`, any number of keyword
arguments may follow. For example

```
    function my_forward_model(sv::RE.RetrievalStateVector; extra_data)
        # Do something with `sv` and `extra_data`
        return true
    end
```

See the online documentation (Types -> Solver Types) for more details.

"""
struct IMAPSolver <: AbstractSolver
    forward_model::Function
    state_vector::AbstractStateVector
    prior_covariance::AbstractMatrix
    max_iterations::Int
    dsigma_scale::AbstractFloat
    dispersions::Dict{<:AbstractSpectralWindow, <:AbstractDispersion}
    indices::Dict{<:AbstractSpectralWindow, <:AbstractVector}
    radiance::Radiance
    jacobians::Dict{<:AbstractStateVectorElement, <:Radiance}
    measured::Dict{<:AbstractDispersion, <:AbstractVector}
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
        jacobians::Dict{<:AbstractStateVectorElement, <:Radiance},
        measured::Dict{<:AbstractDispersion, <:AbstractVector},
        noise::Dict{<:AbstractDispersion, <:AbstractVector}
    )

        # Construct the prior covariance matrix. It's convenient to have it as part of the
        # structure, since it is not changing during the inversion.

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
