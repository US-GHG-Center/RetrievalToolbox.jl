# [Thermal Atmosphere Element Types](@id thermal_types)

The following types can be used to instantiate objects that count as `AbstractAtmosphereElement`, and thus can be used in e.g. `EarthAtmosphere` objects. Their role is to signify to the radiative transfer routines that the surface or the atmosphere itself shall emit thermal radiation, based on Planck's law. If the chosen radiative transfer solver is able, it will automatically set the correct options in the radiative transfer interface. At the moment, only the XRTM RT type is capable of using thermal emission. Users have to make sure to set `source_thermal` in their model options, and only some of XRTM's solvers are able to account for thermal radiation.

For a thermally active atmosphere, the temperature information is taken straight from the `.temperature` profile of the `EarthAtmosphere` object. For a thermal surface, the surface temperature is stated in the object itself. To retrieve the surface temperature via a state vector element, see the [`SurfaceTemperatureSVE`](@ref) state vector element. Retrieving the atmospheric temperature profile is currently possible via the [`TemperatureOffsetSVE`](@ref) state vector element.

```@docs
RE.ThermalAtmosphereIsotropic
```

```@docs
RE.ThermalSurfaceIsotropic
```

## Implementation of Planck's law

```@docs
RE.Planck_radiance(λ::Unitful.Length, T::Unitful.Temperature, rad_unit::Unitful.Units)
```