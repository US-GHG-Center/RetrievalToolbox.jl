"""
$(TYPEDFIELDS)

A type to hold a collection of state vector elements
for use in a retrieval.
"""
struct RetrievalStateVector <: AbstractStateVector

    state_vector_elements::Vector{AbstractStateVectorElement}

end

"""
$(TYPEDFIELDS)

A type to hold a collection of state vector elements
for use in forward model runs (TOA spectrum simulations).
"""
struct ForwardModelStateVector <: AbstractStateVector end


"""
A type to represent a gas profile to be retrieved
"""
mutable struct GasVMRProfileSVE{T} <: AbstractStateVectorElement
    "The atmospheric level affected by this SVE"
    level::Integer
    "The gas object whose VMR profile to be retrieved"
    gas::GasAbsorber
    "Require dimensionless unit (e.g. Unitful.percent, or u\"1\")"
    unit::Unitful.DimensionlessUnits
    "First guess VMR"
    first_guess::T
    "Prior VMR"
    prior_value::T
    "Prior VMR covariance"
    prior_covariance::T
    "Per-iteration VMR values"
    iterations::Vector{T}
end

"""
$(TYPEDFIELDS)

A type to represent a state vector element which scales a gas profile sub-column at the
retrieval grid by its current value. `start_level` must be lower value (higher up in the
atmosphere) than the `end_level`.
"""
mutable struct GasLevelScalingFactorSVE{T} <: AbstractStateVectorElement

    "Start index of the sub-column (inclusive)"
    start_level::Int
    "Stop index of the sub-column (inclusive)"
    end_level::Int
    "The gas object whose volume mixing ratio is scaled"
    gas::GasAbsorber

    "Require dimensionless unit (e.g. Unitful.percent, or u\"1\")"
    unit::Unitful.DimensionlessUnits

    "First guess value"
    first_guess::T
    "Prior value"
    prior_value::T
    "Prior covariance value"
    prior_covariance::T
    "Vector to hold per-iteration values"
    iterations::Vector{T}

    function GasLevelScalingFactorSVE(
        start_level::Int,
        end_level::Int,
        gas,
        unit,
        fg::T,
        prior::T,
        prior_cov::T
        ) where {T}

        @assert start_level <= end_level "Start level must be <= than end level."

        return new{T}(
            start_level,
            end_level,
            gas,
            unit,
            fg,
            prior,
            prior_cov,
            [fg]
        )

    end

end


"""
$(TYPEDFIELDS)

State vector type for Lambertian surface albedo polynomials
"""
mutable struct SurfaceAlbedoPolynomialSVE{
    T1<:AbstractFloat,
    T2<:AbstractSpectralWindow
    } <: AbstractStateVectorElement

    "Spectral window to which this SVE is attached to"
    swin::T2
    "Polynomial coefficient order (0-indexed, order 0 means constant)"
    coefficient_order::Int

    "Depending on the coefficient order \"o\", this unit should be L^{-o}"
    unit::Unitful.Units

    "First guess value"
    first_guess::T1
    "Prior value"
    prior_value::T1
    "Prior covariance value"
    prior_covariance::T1
    "Vector to hold per-iteration values"
    iterations::Vector{T1}

    function SurfaceAlbedoPolynomialSVE(
        swin::T2,
        order::Int,
        supply_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
        fg::T1,
        prior::T1,
        prior_cov::T1
        ) where {T1<:AbstractFloat, T2<:AbstractSpectralWindow}


        # First we check if the supply_unit is compatible with the spectral
        # unit of the spectral window:
        if (supply_unit isa Unitful.LengthUnits) &
            (swin.ww_unit isa Unitful.WavenumberUnits)
            @error "Incompatible spectral units."
            return nothing
        end

        if (supply_unit isa Unitful.WavenumberUnits) &
            (swin.ww_unit isa Unitful.LengthUnits)
            @error "Incompatible spectral units."
            return nothing
        end


        # Dynamically created unit of 1/(Length^order)
        # This is a necessary step in the inner constructor to create a type
        # which will be consistent with the order of the polynomial that this
        # state vector element represents.

        this_unit = supply_unit^-order

        return new{T1, T2}(swin, order, this_unit, fg, prior, prior_cov, [fg])

    end

end

"""
$(TYPEDFIELDS)

State vector type for BRDF amplitude polynomials
"""
mutable struct BRDFPolynomialSVE{
    T1<:AbstractFloat,
    T2<:AbstractSpectralWindow,
    T3<:BRDFKernel
    } <: AbstractStateVectorElement

    "Spectral window to which this SVE is attached to"
    swin::T2
    "BRDF Kernel type on which this SVE acts"
    BRDF_type::Type{T3}
    "Polynomial coefficient order (0-indexed, order 0 means constant)"
    coefficient_order::Int
    "Depending on the coefficient order \"o\", this unit should be L^{-o}"
    unit::Unitful.Units

    "First guess value"
    first_guess::T1
    "Prior value"
    prior_value::T1
    "Prior covariance value"
    prior_covariance::T1
    "Vector to hold per-iteration values"
    iterations::Vector{T1}

    function BRDFPolynomialSVE(
        swin::T2,
        BRDF_type::Type{T3},
        order::Int,
        supply_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
        fg::T1,
        prior::T1,
        prior_cov::T1
        ) where {T1<:AbstractFloat, T2<:AbstractSpectralWindow, T3<:BRDFKernel}


        # First we check if the supply_unit is compatible with the spectral
        # unit of the spectral window:
        if (supply_unit isa Unitful.LengthUnits) &
            (swin.ww_unit isa Unitful.WavenumberUnits)
            @error "Incompatible spectral units."
            return nothing
        end

        if (supply_unit isa Unitful.WavenumberUnits) &
            (swin.ww_unit isa Unitful.LengthUnits)
            @error "Incompatible spectral units."
            return nothing
        end


        # Dynamically created unit of 1/(Length^order)
        # This is a necessary step in the inner constructor to create a type
        # which will be consistent with the order of the polynomial that this
        # state vector element represents.

        this_unit = supply_unit^-order

        return new{T1, T2, T3}(
            swin, BRDF_type, order, this_unit, fg, prior, prior_cov, [fg]
            )

    end

end




"""
$(TYPEDFIELDS)

State vector type for zero level offset polynomials
"""
mutable struct ZeroLevelOffsetPolynomialSVE{
    T1<:AbstractFloat,
    T2<:AbstractSpectralWindow
    } <: AbstractStateVectorElement

    "Spectral window to which this SVE is attached to"
    swin::T2
    "Polynomial coefficient order (0-indexed, order 0 means constant)"
    coefficient_order::Int

    "Depending on the coefficient order \"o\", this unit should be L^{-o}"
    ww_unit::Unitful.Units
    "Full unit"
    unit::Union{Unitful.Units, Real}

    "First guess value"
    first_guess::T1
    "Prior value"
    prior_value::T1
    "Prior covariance value"
    prior_covariance::T1
    "Vector to hold per-iteration values"
    iterations::Vector{T1}

    function ZeroLevelOffsetPolynomialSVE(
        swin::T2,
        order::Int,
        supply_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
        radiance_unit::Union{Unitful.Units, Real},
        fg::T1,
        prior::T1,
        prior_cov::T1
        ) where {T1, T2}

        #=
        Dynamically created unit of 1/(Length^order) This is a necessary step in the inner
        constructor to create a type which will be consistent with the order of the
        polynomial that this state vector element represents.

        NOTE At the moment, for orders > 0, the resulting coefficient unit will be
        "radiance unit / wavelength^order". So when using e.g. "ph/s/m^2/sr/µm" and "µm",
        the unit of the SVE will be "ph/s/m^2/sr/µm^2". Right now, there is no way of
        dictating when Unitful cancels units, so this behavior cannot be changed such that
        the unit will be the "more intuitive" "ph/s/m^2/sr/µm/µm".
        =#

        new_ww_unit = supply_unit^-order # this is a Unitful.Units

        if radiance_unit isa Number
            # Here: radiance_unit * new_ww_unit would be a Unitful.Quantity, so we must
            # extract the unit
            final_ww_unit = unit(radiance_unit * new_ww_unit)
        elseif radiance_unit isa Unitful.Units
            # Here: radiance_unit * new_ww_unit is a Unitful.Units, an we can leave it
            final_ww_unit = radiance_unit * new_ww_unit
        end

        return new{T1, T2}(
            swin, order, new_ww_unit, final_ww_unit,
            fg, prior, prior_cov, [fg]
            )

    end

end



"""
$(TYPEDFIELDS)

State vector type for scaling the solar continuum via a wavelength-dependent
polynomial.
"""
mutable struct SolarScalerPolynomialSVE{T<:AbstractFloat} <: AbstractStateVectorElement

    "Spectral window to which this SVE is attached to"
    swin::AbstractSpectralWindow
    "Polynomial coefficient order (0-indexed, order 0 means constant)"
    coefficient_order::Int
    "Depending on the coefficient order \"o\", this unit should be L^{-o}"
    unit::Unitful.Units
    "First guess value"
    first_guess::T
    "Prior value"
    prior_value::T
    "Prior covariance value"
    prior_covariance::T
    "Vector to hold per-iteration values"
    iterations::Vector{T}

    function SolarScalerPolynomialSVE(
        swin::AbstractSpectralWindow,
        order::Int,
        supply_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
        fg::T,
        prior::T,
        prior_cov::T) where {T}

        this_unit = supply_unit^-order

        return new{T}(swin, order, this_unit, fg, prior, prior_cov, [fg])

    end

end

"""
$(TYPEDFIELDS)

State factor type for dispersion polynomial coefficient retrieval
"""
mutable struct DispersionPolynomialSVE{
    T<:AbstractFloat
    } <: AbstractStateVectorElement

    dispersion::AbstractDispersion
    coefficient_order::Int

    "Require unit of length or wavenumber"
    unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}

    first_guess::T
    prior_value::T
    prior_covariance::T

    iterations::Vector{T}

    function DispersionPolynomialSVE(
        dispersion::AbstractDispersion,
        order::Int,
        unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
        fg::T,
        prior::T,
        prior_cov::T) where {T}

        return new{T}(dispersion, order, unit, fg, prior, prior_cov, [fg])

    end

end

mutable struct ILSStretchPolynomialSVE{T <: AbstractFloat} <: AbstractStateVectorElement

    swin::AbstractSpectralWindow
    coefficient_order::Int
    unit::Unitful.DimensionlessUnits
    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function ILSStretchPolynomialSVE(
        swin::AbstractSpectralWindow,
        order::Int,
        unit::Unitful.DimensionlessUnits,
        fg::T,
        prior::T,
        prior_cov::T
        ) where {T}

        return new{T}(swin, order, unit, fg, prior, prior_cov, [fg])

    end

end

"""
$(TYPEDFIELDS)

State factor type for surface pressure retrieval
"""
mutable struct SurfacePressureSVE{T <: AbstractFloat} <: AbstractStateVectorElement

    "Require pressure units (e.g. Pa, Torr, hPa, ..)"
    unit::Unitful.PressureUnits

    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function SurfacePressureSVE(
        unit::Unitful.PressureUnits,
        fg::T,
        prior::T,
        prior_cov::T,
        ) where {T}

        return new{T}(unit, fg, prior, prior_cov, [fg])

    end

end

"""
$(TYPEDFIELDS)

State factor type for temperature profile (constant) offset
"""
mutable struct TemperatureOffsetSVE{T <: AbstractFloat} <: AbstractStateVectorElement

    unit::Unitful.Units{U, Unitful.𝚯, nothing} where U

    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function TemperatureOffsetSVE(
        unit::Unitful.Units{U, Unitful.𝚯, nothing} where U,
        fg::T,
        prior::T,
        prior_cov::T
    ) where {T}

        return new{T}(unit, fg, prior, prior_cov, [fg])

    end

end

"""
$(TYPEDFIELDS)

State factor type for SIF radiance. This controls the `radiance_at_reference` field of the
`SIF` object that is referenced here.
"""
mutable struct SIFRadianceSVE{T <: AbstractFloat} <: AbstractStateVectorElement

    SIF::AbstractSIFRadiance
    unit::Unitful.Units
    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function SIFRadianceSVE(
        SIF::AbstractSIFRadiance,
        fg::T,
        prior::T,
        prior_cov::T
    ) where {T}

        return new{T}(SIF, SIF.radiance_unit, fg, prior, prior_cov, [fg])

    end

end


"""
NOTE! This is not implemented yet anywhere!
"""
mutable struct SolarSpectralShiftSVE{T} <: AbstractStateVectorElement

    spectral_window::AbstractSpectralWindow
    unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function SolarSpectralShiftSVE(
        swin::AbstractSpectralWindow,
        unit,
        fg::T,
        prior::T,
        prior_cov::T,
    ) where {T}

        return new{T}(swin, unit, fg, prior, prior_cov, [fg])

    end

end

mutable struct AerosolOpticalDepthSVE{T} <: AbstractStateVectorElement

    aerosol::AbstractAerosolType
    log::Bool
    unit::Unitful.DimensionlessUnits
    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function AerosolOpticalDepthSVE(
        aerosol,
        log,
        unit,
        fg::T,
        prior::T,
        prior_cov::T
    ) where {T}

        return new{T}(aerosol, log, unit, fg, prior, prior_cov, [fg])

    end

end

mutable struct AerosolHeightSVE{T} <: AbstractStateVectorElement

    aerosol::AbstractAerosolType
    log::Bool
    unit::Unitful.DimensionlessUnits
    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function AerosolHeightSVE(
        aerosol,
        log,
        unit,
        fg::T,
        prior::T,
        prior_cov::T
    ) where {T}

        return new{T}(aerosol, log, unit, fg, prior, prior_cov, [fg])

    end

end

mutable struct AerosolWidthSVE{T} <: AbstractStateVectorElement

    aerosol::AbstractAerosolType
    log::Bool
    unit::Unitful.DimensionlessUnits
    first_guess::T
    prior_value::T
    prior_covariance::T
    iterations::Vector{T}

    function AerosolWidthSVE(
        aerosol,
        log,
        unit,
        fg::T,
        prior::T,
        prior_cov::T
    ) where {T}

        return new{T}(aerosol, log, unit, fg, prior, prior_cov, [fg])

    end

end


#=

    "Accessor"-type functions

=#


get_unit(SVE::AbstractStateVectorElement) = SVE.unit
get_unit(SV::AbstractStateVector) = map(
    get_unit, SV.state_vector_elements)

get_first_guess(SVE::AbstractStateVectorElement) = SVE.first_guess
get_first_guess_with_unit(SVE::AbstractStateVectorElement) =
    get_first_guess(SVE) * SVE.unit
get_current_value(SVE::AbstractStateVectorElement) = SVE.iterations[end]
get_current_value_with_unit(SVE::AbstractStateVectorElement) =
    get_current_value(SVE) * SVE.unit
get_prior_value(SVE::AbstractStateVectorElement) = SVE.prior_value
get_prior_value_with_unit(SVE::AbstractStateVectorElement) =
    get_prior_value(SVE) * SVE.unit
get_prior_covariance(SVE::AbstractStateVectorElement) = SVE.prior_covariance
get_prior_covariance_with_unit(SVE::AbstractStateVectorElement) =
    get_prior_covariance(SVE) * SVE.unit^2
get_prior_uncertainty(SVE::AbstractStateVectorElement) = sqrt(get_prior_covariance(SVE))


#=
From here on, we dynamically create appropriate functions to
return the above quantities from the full state vector, rather
than just a state vector element
=#

for fun in (
    :get_first_guess,
    :get_current_value,
    :get_prior_value,
    :get_prior_covariance,
    :get_prior_uncertainty
)

    for _u in ["_with_unit", ""]

        fstr = "$(String(fun))$(_u)(SV::AbstractStateVector) = " *
               "map($(String(fun))$(_u), SV.state_vector_elements)"

        eval(Meta.parse(fstr))

    end

end