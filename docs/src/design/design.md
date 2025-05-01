# Toolkit design philosophy

## Lack of exported functions (namespace)



## Reliance on Julia types

## Buffers for performance

## Executing the forward model will mutate objects

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
swin.wavelength;
show(swin.wavelength')
```

or even the Unicode symbol λ:
```@repl swin
swin.λ;
show(swin.λ')
```

Note that `swin.λ` or `swin.wavelength` do not perform a calculation or conversion. For the spectral window type `SpectralWindow`, an overloaded `getproperty` function was dynamically created during startup such that `getproperty(SpectralWindow, :λ)` returns the `ww` field, and similarly for wavenumber units.

Since types with some spectral dimension must also have a corresponding unit field, `ww_unit`, the `getproperty` function is able to check whether the requested spectral unit is appropriate. Trying to access `swin.wavenumber` or `swin.ν` will fail:

```@repl swin
swin.wavenumber # or swin.ν
```

Benchmarks have shown there is no significant performance drawback due to the use of the new `getproperty` functions.

When writing new functions that use objects of any types with a spectral unit, users can interrogate the `ww_unit` field to control the behavior of calculations. Inside RetrievalToolbox, this is done, for example, for the Doppler shift calculations, which use different expressions depending on whether they are performed in wavelength or wavenumber space.

```@example swin
function something_new(swin::RE.SpectralWindow)

    if swin.ww_unit isa Unitful.LengthUnit
        # Do calculations in wavelength space
    elseif swin.ww_unit isa Unitful.WavenumberUnit
        # Do calculations in wavenumber space
    end

end;
```

When users write functions, they should be mindful of when they use `.ww` compared to `.wavelength` or `.wavenumber`. For many calculations, it makes no difference whether the spectral unit is wavelength or wavenumber, and thus writing `.ww` is legible and reasonable. If some function only makes sense in one spectral unit, but not in the other, then `.wavenumber` or `.wavelength` can be used. Note that in this case, an error will be thrown if an object with the wrong spectral unit is passed into this function, and the spectral unit will be accessed via `.wavenumber` or `.wavelength`. This might be the desired behavior - calculations should throw an error if invoked in the wrong spectral unit space.

When users write top-level retrieval scripts for some specific scenario, it is usually most obvious to write the specific spectral unit, as that does not change and it also becomes clear that only that specific spectral unit is considered. I.e., a retrieval script that launches OCO-2 retrievals should access the spectral unit with `.wavelength`.


## Build your own algorithm!

### Custom forward model

### Lack of an instrument type

### Custom ingestion of needed inputs

## Extend the toolkit with your own types and functions
