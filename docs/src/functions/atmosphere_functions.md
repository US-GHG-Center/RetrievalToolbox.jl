# Atmosphere functions

Functions below are designed to operate on atmosphere object, and for the time being
mostly on `EarthAtmosphere` ones.

## Updates and roll-backs

Atmosphere objects, just as `EarthAtmosphere`, are generally mutated during the course of an inverse retrieval problem. There is no "history" per se that would allow to unambiguously re-create the atmospheric state at some earlier iteration. Due to the flexibility of the toolset, it is always possible to make changes (intentional or not) to the atmospheric state that would be difficult to revert.

The mechanism that is used in RetrievalToolbox relies on user discipline. Let `atm` be an `EarthAtmosphere` that has been initialized and adjusted to the user's wishes, and let `SV` be a `RetrievalStateVector` that contains state vector elements that imply some modification of the contents of `atm`.

For some parts of the atmosphere, their state can be easily and directly inferred from a state vector element. As an example we can look at `SurfacePressureSVE`: the lowest pressure level is simply set to the *current* value of that state vector element via the `atmosphere_element_statevector_update!` function. One can think of the surface pressure of the model atmosphere being quite rigidly controlled, it will simply take on whatever value `SurfacePressureSVE` currently has. This allows us to also recreate the surface pressure that was present at an earlier iteration. We could manipulate the contents of the `.iterations` field of the `SurfacePressureSVE` object, i.e. remove entries up to the requested iteration number. When we call `atmosphere_element_statevector_update!` again, the atmosphere will revert back to that earlier state.

The above holds true for many state vector elements that control the atmospheric state, like the `GasVMRProfileSVE` which defines the volume mixing ratio profile of a certain gas in that model atmosphere, or `SIFRadianceSVE` which defines the magnitude of surface-leaving fluorescence radiance emission.

There are notable exceptions, for which this concept does not work - at least not reliably. `GasLevelScalingFactorSVE`, for example, scales the *initial* gas profile by some factor. Similarly, `TemperatureOffsetSVE` adds some value to the *initial* temperature profile. While it is mathematically straightforward to reconstruct the these profiles from the implicit history of iterations that is found in the `.iterations` vectors of these state vector elements, it has proven to be not fully reliable in practice.

Instead, RetrievalToolbox expects users to to call a function `atmosphere_element_statevector_rollback!` towards the end of their forward model implementation, which re-sets particular atmospheric components to their initial state. Therefore, when a new forward model evaluation is performed, the call to `atmosphere_element_statevector_update!` will modify those components correctly to the desired state as given by the state vector elements.

!!! warning
    **Important!** There is no mechanism to check if the *update* and *rollback* functions have been called in a forward model. It is fully up to users to ensure that those calls are made in the forward model in the appropriate places.


To summarize, users should make a call to `atmosphere_element_statevector_update!` towards the beginning of their forward model, and then another call to `atmosphere_element_statevector_rollback!` towards the end.

```julia
# My forward model

# [...]
for atm_e in atm.atm_elements
    RE.atmosphere_element_statevector_update!(atm_e, SV)
end

# [...]
# calculate radiances, Jacobians etc.
# [...]

for atm_e in atm.atm_elements
    RE.atmosphere_element_statevector_rollback!(atm_e, SV)
end

# [...]
# forward model end
```

!!! tip
    Users also must ensure that `atmosphere_element_statevector_update!` is called *after* iterations are complete if they want to e.g. calculate some atmospheric quantity corresponding to the final state inferred by the retrieval.


Note that the code snipped shows explicit loops over all atmosphere elements `atm_e` of some atmosphere `atm`.


## List of atmosphere functions
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["atmosphere.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```