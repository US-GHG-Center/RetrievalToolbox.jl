# Toolkit design philosophy

## Namespace

The RetrievalToolbox module exports many functions, types and variables - some of those might share names with functions from other modules or your own user code. While not strictly necessary, we generally recommend to load the module and then declare an alias to then call the module functions through it.

```julia
using RetrievalToolbox
const RE = RetrievalToolbox

RE.some_function()

# this, however, works too
some_function()
```

## Reliance on Julia types

One of Julia's most prominent features is the flexibility that comes with its [rich type system](https://docs.julialang.org/en/v1/manual/types/). While not object-oriented in the sense that C++ is, objects, their relationship to functions, and how they act on specific objects is a major design component of Julia and thus RetrievalToolbox.

RetrievalToolbox defines a number of abstract types and then lots of composite types, which are akin to `struct` types in C. They usually sit below abstract types in the type hierarchy and represent some granular concept that is useful for trace gas retrievals.

For example, the `GaussAerosol` type belongs to the `AbstractAerosolType` which itself belongs to `AbstractAtmosphereElement`, and describes some aerosol whose vertical distribution in the model atmosphere is described by a Gaussian. The full type hierarchy here is `GaussAerosol ⊂ AbstractAerosolType ⊂ AbstractAtmosphereElement`.
When users create a model atmosphere, they must add some `AbstractAtmosphereElement` to the list of atmosphere elements. At that stage, we do not care what _specific_ object that might be, as long as it satisfies the requirement that it has to be a type that is a subtype of `AbstractAtmosphereElement`.

Now how does Julia's type system and the **multiple dispatch** paradigm help us here. First, it allows for some convenience. We can use some list of atmospheric elements that we want represented in our model atmosphere: aerosols and Rayleigh scattering. Those two are quite different in practical terms, even though they might act in similar ways on our various calculations. So if we have some list of `AbstractAtmosphereElement` objects

```julia
atm_list = [aerosol1, aerosol2, RayleighScattering()]
```

we would ideally want to perform some action with each of those elements, such as calculating their contribution to the optical depth profiles of our model atmosphere. Some naïve way of doing that would be to iterate through each element and perform the appropriate action:

```julia
for atm in atm_list
    if is_an_aerosol(atm)
        tau = calculate_gauss_aerosol_tau(atm)
    end

    if is_Rayleigh_scattering(atm)
        tau = calculate_Rayleigh_tau(atm)
    end

    # .. do something with tau
end
```

This is easily done with Julia, but makes this top-level iteration not very elegant. One can imagine that if we implement several aerosol distribution types, this loop will grow since we want to call the dedicated function to perform the wanted operation. The preferred way in Julia is making use of multiple dispatch: we decide on a name for a function that shall perform the equivalent task for different types. Let us call this function `calculate_tau`:

```julia
for atm in atm_list
    tau = calculate_tau(atm)
    # .. do something with tau
end
```

and as long as there is a function `calculate_tau` which implements the calculate for the requested type, above code will execute and keep this top-level loop nice and tidy.

While the above example allows for some convenience, the strength of multiple dispatch also lies in how users can expand code without having to change code deep within a module they use. For example, let us imagine a new aerosol type that a user wants to integrate in their retrieval application: `MyAerosolType`. A new list of atmosphere elements would be created, like so

```julia
my_new_aerosol = MyAerosolType()
atm_list = [aerosol1, aerosol2, RayleighScattering(), my_new_aerosol()]
```

Now in this imagined scenario, the user can now create their own `calculate_tau` function, which would implement the specific routine that computes the optical depth profiles for their new aerosol type:

```julia
function calculate_tau(a::MyAerosolType)
    # Do many calculations here...
    # ...
end
```

When done correctly, the new function will be invoked when the loop above runs (`for atm in atm_list ...`)  since the Julia compiler will now look for a `calculate_tau` function that can act on an object of type `MyAerosolType`.

To summarize: RetrievalToolbox makes extensive use of Julia's type system such that many functions in the module do not act on primitive types (like numbers or strings), but on custom composite types.

!!! info
    The Julia type system allows for flexibility in creating program code that looks the same for a variety of object types. Further, users can more easily extend existing code with their own types without necessarily having to change the underlying routines, but by providing their own user code.


## Buffers for performance

Julia is a garbage-collected language (more detailed info [here](https://docs.julialang.org/en/v1/devdocs/gc/)), meaning that users do not have explicit control over how and when objects are de-allocated from memory. When users write functions that allocate (by creating vectors and arrays, for example), the memory is not immediately freed when the function is completed. The garbage collector (GC) is triggered at some point when the memory usage reaches some level. The GC then traces through the objects in memory and removes those that are no longer used.

The big advantage of GC-based languages is of course that manual memory management is no longer needed, and users do not have to keep track of correct allocation and de-allocation of objects, and memory safety issues are also less common. The major downside is that without manual memory management, users could (un)willingly write code that very inefficiently allocates a lot of memory. When those allocations happen in certain places (loops mostly), memory will fill up quickly and trigger GC sweeps very often. The paradigm within Julia is usually: allocate arrays and vectors beforehand, and perform calculations on these pre-allocated objects.

Early versions of RetrievalToolbox did not make use of much pre-allocation, and most calls to functions would create new objects. This has proven to be not a feasible solution. While convenient for top-level scripting, the many allocations needed made it impossible to run faster retrievals where the forward model run was much less than a second (e.g. physics-based SIF retrievals or other non-scattering applications).

Thus, we make use of pre-allocated objects, which we call buffers here.

## Executing the forward model will mutate some objects

Kinda bad for e.g. gas scale factors that mean factors of some initial atmospheric state.

## Wrapper functions and specific dispatch

## Considering quantities with physical units

Explain types and dedicated type fields for units. Pay attention to supplying quantities with the right units!

## Wavelengths and wavenumbers

RetrievalToolbox supports two fundamental spectral unit types: wavelength and wavenumber. Users might want to build an algorithm pipeline that is specific to some instrument, which natively produces spectra in either wavelength or wavenumber units. In order to make the native spectral unit be visible as such, RetrievalToolbox provides dynamic accessors which allow users to reference the spectral unit of objects using their natural wording or symbol.

Rather than writing duplicate types and functions, RetrievalToolbox employs _magic accessor_ methods. Any quantity that represents a spectral unit, is typed `ww` or, for example, `ww_unit` or `ww_reference`. The two-letter combination `ww` is thus reserved in the RetrievalToolbox codebase, and no type fields should contain this combination of letters.

The `ww` should be considered a placeholder, which represents either a wavelength- or a wavenumber-related quantity. Any type that contains a field or quantity `ww` must also contain a field named `ww_unit`. Other fields are optional.

When the RetrievalToolbox module is imported, all types inside the RetrievalToolbox namespace are scanned for type fields that contain the substring `ww`. For each type that contains such a field, a new accessor function is dynamically created, which allows users to access the spectral type fields with the appropriate symbol. Illustrative examples follow.

Creating a spectral window object from 1.49 µm through 1.55 µm with 10 nm spacing could, for example, look like this (with loaded `Unitful`):

```@repl swin; continued = true
using RetrievalToolbox # hide
const RE = RetrievalToolbox # hide
using Unitful # hide

swin = RE.SpectralWindow(
    "test",
    1.49, # Lower limit
    1.55, # Upper limit
    collect(1.48:0.01:1.56), # Create the grid with spacing
    Unitful.µm, # Designate microns as unit of choice
    1.50 # Set some reference wavelength
    )
```

As can be seen in the type definition, the spectral grid can be accessed via `swin.ww_grid`.

```@repl swin
swin.ww_grid;
show(swin.ww_grid')
```

Now the magic accessor allows users to access the same field using the more "natural" wavelength term

```@repl swin
swin.wavelength_grid;
show(swin.wavelength_grid')
```

or even the Unicode symbol λ:
```@repl swin
swin.λ_grid;
show(swin.λ_grid')
```

Note that `swin.λ_grid` or `swin.wavelength_grid` do not perform a calculation or conversion. For the spectral window type `SpectralWindow`, an overloaded `getproperty` function was dynamically created during startup such that `getproperty(SpectralWindow, :λ)` returns the `ww_grid` field, and similarly for wavenumber units.

Since types with some spectral dimension must also have a corresponding unit field, `ww_unit`, the `getproperty` function is able to check whether the requested spectral unit is appropriate. Trying to access `swin.wavenumber_grid` or `swin.ν_grid` will fail:

```@repl swin
swin.wavenumber_grid # or swin.ν
```

Benchmarks have shown there is no significant performance drawback due to the use of the new `getproperty` functions.

When writing new functions that use objects of any types with a spectral unit, users can interrogate the `ww_unit` field to control the behavior of calculations. Inside RetrievalToolbox, this is done, for example, for the Doppler shift calculations, which use different expressions depending on whether they are performed in wavelength or wavenumber space.

```julia
function something_new(swin::RE.SpectralWindow)

    if swin.ww_unit isa Unitful.LengthUnit
        # Do calculations in wavelength space
    elseif swin.ww_unit isa Unitful.WavenumberUnit
        # Do calculations in wavenumber space
    end

end
```

When users write functions, they should be mindful of when they use `.ww_grid` compared to `.wavelength_grid` or `.wavenumber_grid`. For many calculations, it makes no difference whether the spectral unit is wavelength or wavenumber, and thus writing `.ww` is legible and reasonable. If some function only makes sense in one spectral unit, but not in the other, then `.wavenumber_grid` or `.wavelength_grid` can be used. Note that in this case, an error will be thrown if an object with the wrong spectral unit is passed into this function, and the spectral unit will be accessed via `.wavenumber_grid` or `.wavelength_grid`. This might be the desired behavior - calculations should throw an error if invoked in the wrong spectral unit space.

When users write top-level retrieval scripts for some specific scenario, it is usually most obvious to write the specific spectral unit, as that does not change and it also becomes clear that only that specific spectral unit is considered. I.e., a retrieval script that launches OCO-2 retrievals should access the spectral unit with `.wavelength_grid`.


## Build your own algorithm!

### Custom forward model

### Lack of an instrument type

### Custom ingestion of needed inputs

## Extend the toolkit with your own types and functions
