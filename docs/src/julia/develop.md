# [How to develop or extend RetrievalToolbox](@id develop_RetrievalToolbox)

Thanks to Julia's approach to user-defined types and functions that act on or use those types, it is very straightforward to extend the capabilities of RetrievalToolbox to users' own needs, without having to alter the RetrievalToolbox library itself.

## Extending RetrievalToolbox

Users can extend RetrievalToolbox for use in their own projects by adding or overwriting its existing functions. If wanted, these new functions can seamlessly propagate into the overall flow of functions *within* RetrievalToolbox. The following section will provide some examples that can be used as guidance.

### Example 1: A user function that involves a RetrievalToolbox type

In this first example, we will simply create a new function to extend the functionality of RetrievalToolbox. This new function will make use of existing functions and types that RetrievalToolbox provides.

Let the following be the task at hand:

!!! tip "Task"
    We want to write a function which calculates the total, isotropic solar-induced fluorescence (SIF) emission into the half hemisphere (towards the sky). The function will take a `SIFRadiance` type as argument, and simply produce the integral over the spectral dimensions and the positive half space.


We start off by creating a `SIFRadiance` object and just inspect the spectrally resolved isotropic radiance that it represents. Using the function `get_SIF_radiance`, we can simply obtain that emitted radiance for some given wavelength (or wavenumber).

```@example sif
using RetrievalToolbox; const RE = RetrievalToolbox
using Unitful
using Plots; default(titlefontsize=10, labelfontsize=8)

# Create the SIF object with 1e-7 W/m2/sr/µm at a reference wavelength of 740nm.
sif = RE.SIFRadiance(1.0e-7, u"W/m^2/sr/μm", 740.0, u"nm");

# Make a plot for sake of visualization
wl = collect(650.0u"nm":1.0u"nm":870u"nm")
Plots.plot(
    wl,
    RE.get_SIF_radiance.(Ref(sif), wl),
    label=nothing, size=(400, 200)
)
xlabel!("Wavelength")
ylabel!("SIF radiance\n($(sif.radiance_unit))")
Plots.savefig("develop_fig1.png"); nothing # hide
```
![sif plot](develop_fig1.png)


To get the spectrally integrated value, we only need to sample this emitted radiance waveform in spectral space and integrate with some appropriate technique (trapezoidal integration is sufficient here). Further, we then multiply the spectrally integrated value with ``2\pi`` to take care of the integration of the half-space, which is equal to ``2\pi`` steradians. Also note how we use various units to make sure that we get the correct units for the final result in Watts per square meter.

The function only takes one argument, `sif`, which must be of the sub-type `AbstractSIFRadiance`, which RetrievalToolbox exports.

```@example sif
function SIF_integral(sif::T) where T<:RE.AbstractSIFRadiance

    # Spectrally integrated value
    int_value = 0. * u"W/m^2/sr"
    for i in 2:length(sif.ww_waveform)

        # Simply sample the waveform data underlying
        # the `sif` object..
        ww1 = sif.ww_waveform[i]
        ww2 = sif.ww_waveform[i-1]

        v1 = RE.get_SIF_radiance(sif, ww1) * sif.radiance_unit
        v2 = RE.get_SIF_radiance(sif, ww2) * sif.radiance_unit

        int_value += 0.5 * (v1 + v2) * (ww1 - ww2)
    end

    # Integrate over half-space (2π sr),
    # we have to "manually" set the units to W/m^2 since
    # `Unitful` likes to turn W into kg/m/s^2.
    return (int_value * (2*pi*u"sr") |> u"W/m^2")
end
```

And we finally call the function to test it:

```@example sif
# Amount emitted towards the sky per m^2 of surface area for the
# `sif` object defined above.
SIF_integral(sif)
```

Users can write a function like the above in a script, for example, and simply use it as needed in that particular `.jl` file. If users want something more lasting, since this might be a function they would use over and over again, they could write a new module and use that module in their applications:

```julia
module MyNewModule

    # import other module
    using RetrievalToolbox
    const RE = RetrievalToolbox

    export SIF_integral # Make this function available to users of `MyNewModule`

    # Add SIF function
    function SIF_integral(sif::T) where T<:RE.AbstractSIFRadiance
        # [...] relevant code
    end

end
```

In the example above, our function would be accessible via `MyNewModule.SIF_radiance`.

To summarize, we can easily add new functions that make use of existing RetrievalToolbox functions and thus extend the overall capability of the software library to meet our own specific needs. In this example, we do not need to change anything within RetrievalToolbox itself, we simple make use of the existing functions. Further, we can create new modules that provide our new function if we need to use them for many projects.


### Example 2: A new user type that will become part of the existing RetrievalToolbox type hierarchy

A step up in complexity is the introduction of a new user type. For example, if some form of atmospheric constituent or atmospheric element is needed that is currently not supported in RetrievalToolbox.

In RetrievalToolbox, we have an existing type hierarchy that starts with `AbstractAtmosphereElement`, and then moves on to another abstract type. For example, the SIF object we have used in Example 1 is of type `SIFRadiance`, which itself is a type that belongs to the abstract type `AbstractSIFRadiance` that is an `AbstractAtmosphereElement`.

Below is a visual representation of the current type hierarchy for all `AbstractAtmosphereElements`.

```@example sif
using InteractiveUtils
using AbstractTrees
AbstractTrees.children(x::Type) = subtypes(x)
print_tree(RE.AbstractAtmosphereElement)
```

!!! tip "Task"
    A scenario that a user might face is the following. The currently implemented method to calculate the emitted SIF radiance is not exactly what some user wants, but instead prefers their own formulation. The task here is therefore to create a new, user-defined object that represents our own SIF, and then implement the appropriate functions that compute the radiance.

For the sake of simplicity, let us just assume our new SIF produces a constant radiance of ``1\mathrm{W}/\mathrm{m}^2/\mathrm{sr}/\mathrm{µm}`` between 700 nm and 750 nm. We start by defining our new type, which is similar to the `SIFRadiance` type found in `src/types/sif_types.jl`, however omits several fields: `radiance_at_reference` and `ww_reference` (no need to impose some spectrally dependent scaling), and `waveform` and `itp` fields (we have no need for the interpolation object that captures the shape of the SIF waveform). Also, we only allow wavelengths here (again, for simplicity). Our new SIF type shall be called `MySIFRadiance`:

```@example sif
mutable struct MySIFRadiance <: RE.AbstractSIFRadiance

    radiance_unit::Unitful.Units
    ww_unit :: Unitful.LengthUnits
    ww_waveform :: Vector

end
```


Along with the type definition, it is also important to define an external constructor function:

```@example sif
function MySIFRadiance(
    radiance_unit::Unitful.Units,
    ww_unit::Unitful.LengthUnits
)

    # Produce "waveform" grid
    ww_waveform = [700.0u"nm", 750.0u"nm"] .|> ww_unit

    return MySIFRadiance(
        radiance_unit,
        ww_unit,
        ww_waveform
    )

end
```

With both these defined, we can create a new SIF object, which we shall call `new_sif`:

```@example sif
new_sif = MySIFRadiance(u"W/m^2/sr/μm", u"nm");
```

Now, RetrievalToolbox provides a `get_SIF_radiance` function that we used earlier evaluate the SIF radiance at some wavelength or wavenumber. Let's try to use it here to get the SIF radiance (of `new_sif`, which is an object of our new type `MySIFRadiance`) at 725 nm:

```@repl sif
RE.get_SIF_radiance(new_sif, 725.0u"nm")
```

The above call fails: there is no `get_SIF_radiance` function that accepts a SIF object of type `MySIFRadiance`! This is very much intended behavior, since the `get_SIF_radiance` function that RetrievalToolbox provides only makes sense for `SIFRadiance` objects. We have to write our own! While generally we are free to choose the arguments for our own implementation, it is considered **good practice** to **keep the same arguments**. That way, our new user-defined function have a good chance of seamlessly integrating into function inside of RetrievalToolbox.

```@example sif
import RetrievalToolbox: get_SIF_radiance

function RetrievalToolbox.get_SIF_radiance(
    sif::MySIFRadiance, # Note that this is now `MySIFRadiance`
    ww::Unitful.Length
)

    if (ww < sif.ww_waveform[1]) | (ww > sif.ww_waveform[end])
        return 0.
    else
        return 1.0
    end

end
```

Now that we have created a function that accepts our new SIF type `MySIFRadiance`, a call to `get_SIF_radiance` will use the correct method!

```@example sif
RE.get_SIF_radiance(new_sif, 725.0u"nm")
```

Further, any uses of `SIF_integral` that are invoked with objects of the new SIF type `MySIFRadiance` will themselves make calls to the correct `get_SIF_radiance`. **This is an important part of this example!** We did not have to make any adjustments to `SIF_integral` at all, since we kept the same function signature in our new implementation of `get_SIF_radiance`. Julia recognizes, that we call `SIF_integral` with a new argument type (that being `MySIFRadiance`) and re-compiles the function.

```@example sif
SIF_integral(new_sif)
```

Here is a summary of the important steps when creating a new type that fits into the existing type hierarchy.

!!! note "Summary"
    1. Investigate the RetrievalToolbox type hierarchy to understand where your new type fits in best.
    2. Write the new type definition, making sure that you understand which type fields are potentially needed by other functions, or if they can safely omitted if they serve no purpose for your new type.
    3. Write all needed functions to implement the desired behavior for your new type
        a. This requires knowing which functions may use your new type. The Julia function `methodswith` can be useful in finding that out. E.g. `methodswith(RE.SIFRadiance, RetrievalToolbox)` will list all functions that have at least one `SIFRadiance` function argument.
    4. Run tests!


