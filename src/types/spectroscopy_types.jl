"""
$(TYPEDEF)
Empty spectrscopy object.
"""
struct NoSpectroscopy <: AbstractSpectroscopy end


"""
$(TYPEDEF)

Type to hold a full 4D (OCO-type) ABSCO spectroscopy table, with absorption coefficients
as a function of a broadener gas (H2O usually), pressure, temperature, and wavelength.

$(TYPEDFIELDS)

Note that within a typical ABSCO data file, the wavenumber array is (generally) an H5
double type, and all other arrays are likely 32 bit floats. However, we want to be as
general as possible, so we use four different types T1 through T4 to allow for any mixture
of number types for the axes.
"""
struct ABSCOSpectroscopy4D{
    T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real,
    A <: AbstractArray{T4, 4}
    } <: ABSCOSpectroscopy

    "Filename/path of the underlying ABSCO data"
    file_name::String
    "Name of the gas whose cross sections this object represents"
    gas_name::String
    "User-defined factor that was used to scale the entire table"
    scale_factor::T1
    "Spectral axis (wavelength or wavenumber)"
    ww::Vector{T2}
    "Unit belonging to spectral axis (length or wavenumber unit)"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "Temperature axis (2-dimensional)"
    temperatures::Array{T3, 2}
    "Temperature unit"
    temperatures_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U
    "Pressure axis (1-dimensional)"
    pressures::Vector{T3}
    "Pressure unit"
    pressures_unit::Unitful.PressureUnits
    "Broadener volume mixing ratio (unit 1)"
    broadener_vmrs::Vector{T3}
    "Cross section array, order: λ (or ν), H2O, T, p"
    cross_section::A
    "Cross section unit (usually cm^2/molecule)"
    cross_section_unit::Unitful.AreaUnits

end

"""
$(TYPEDEF)

Type to hold a full 3D (OCO-type) ABSCO spectroscopy table, with absorption coefficients
as a function of pressure, temperature, and wavelength.

$(TYPEDFIELDS)

Note that within a typical ABSCO data file, the wavenumber array is (generally) an H5
double type, and all other arrays are 32 bit floats, hence why there are two types T1,T2
for for this structure.
"""
struct ABSCOSpectroscopy3D{
    T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real,
    A <: AbstractArray{T4, 3}
    } <: ABSCOSpectroscopy

    "Filename/path of the underlying ABSCO data"
    file_name::String
    "Name of the gas whose cross sections this object represents"
    gas_name::String
    "User-defined factor that was used to scale the entire table"
    scale_factor::T1
    "Spectral axis (wavelength or wavenumber)"
    ww::Vector{T2}
    "Unit belonging to spectral axis (length or wavenumber unit)"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "Temperature axis (2-dimensional)"
    temperatures::Array{T3, 2}
    "Temperature unit"
    temperatures_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U
    "Pressure axis (1-dimensional)"
    pressures::Vector{T3}
    "Pressure unit"
    pressures_unit::Unitful.PressureUnits
    "Cross section array, order: λ (or ν), T, p"
    cross_section::A
    "Cross section unit (usually cm^2/molecule)"
    cross_section_unit::Unitful.AreaUnits

end

"""
$(TYPEDEF)

Type to hold a full 4D (AER-type) ABSCO spectroscopy table, with absorption coefficients
as a function of a broadener gas (H2O usually), pressure, temperature, and wavelength.
This type is meant to hold spectroscopy data generated via
https://github.com/ReFRACtor/ABSCO.

$(TYPEDFIELDS)

Note! While looking quite similar to an OCO-type ABSCO table, there are some major
differences, such as the dimension order in the cross section array. This type is meant
to be used with outputs of the ABSCO program available at
https://github.com/ReFRACtor/ABSCO, without any further modification! Most notably, this
type of spectroscopy is the result of a calculation that takes place between pressure
levels, but is evaluated at specific temperatures and broadener VMRs. As such, we can
user interpolation in the T and H2O dimension, but do not interpolate int the pressure
dimension.
"""
struct ABSCOAERSpectroscopy4D{
    T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real,
    A <: AbstractArray{T4, 4}
} <: ABSCOSpectroscopy

    "Filename/path of the underlying ABSCO data"
    file_name::String
    "Name of the gas whose cross sections this object represents"
    gas_name::String
    "User-defined factor that was used to scale the entire table"
    scale_factor::T1
    "Spectral axis (wavelength or wavenumber)"
    ww::Vector{T2}
    "Unit belonging to spectral axis (length or wavenumber unit)"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "Temperature array"
    temperatures::Array{T3, 2}
    "Temperature unit"
    temperatures_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U
    "Pressures"
    pressures::Vector{T3}
    "Pressure unit"
    pressures_unit::Unitful.PressureUnits
    "Broadener volume mixing ratio (unit 1)"
    broadener_vmrs::Vector{T3}
    "Cross section array, order: H2O, p, T, λ (or ν)"
    cross_section::A
    "Cross section unit (usually cm^2/molecule)"
    cross_section_unit::Unitful.AreaUnits

end


"""
$(TYPEDEF)

Type to hold a full 4D (AER-type) ABSCO spectroscopy table, with absorption coefficients
as a function of a broadener gas (H2O usually), pressure, temperature, and wavelength.
This type is meant to hold spectroscopy data generated via
https://github.com/ReFRACtor/ABSCO.

$(TYPEDFIELDS)

Note! While looking quite similar to an OCO-type ABSCO table, there are some major
differences, such as the dimension order in the cross section array. This type is meant
to be used with outputs of the ABSCO program available at
https://github.com/ReFRACtor/ABSCO, without any further modification! Most notably, this
type of spectroscopy is the result of a calculation that takes place between pressure
levels, but is evaluated at specific temperatures and broadener VMRs. As such, we can
user interpolation in the T and H2O dimension, but do not interpolate int the pressure
dimension.
"""
struct ABSCOAERSpectroscopy3D{
    T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real,
    A <: AbstractArray{T4, 3}
} <: ABSCOSpectroscopy

    "Filename/path of the underlying ABSCO data"
    file_name::String
    "Name of the gas whose cross sections this object represents"
    gas_name::String
    "User-defined factor that was used to scale the entire table"
    scale_factor::T1
    "Spectral axis (wavelength or wavenumber)"
    ww::Vector{T2}
    "Unit belonging to spectral axis (length or wavenumber unit)"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "Temperature array"
    temperatures::Array{T3, 2}
    "Temperature unit"
    temperatures_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U
    "Pressures"
    pressures::Vector{T3}
    "Pressure unit"
    pressures_unit::Unitful.PressureUnits
    "Cross section array, order: p, T, λ (or ν)"
    cross_section::A
    "Cross section unit (usually cm^2/molecule)"
    cross_section_unit::Unitful.AreaUnits

end


# Not implemented yet
struct HITRANSpectroscopy{T1,T2,T3} <: AbstractSpectroscopy

    gas_name::String
    scale_factor::T1

    ww::Vector{T2}
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}

    temperatures::Array{T1, 2}
    temperatures_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U

    pressures::Vector{T1}
    pressures_unit::Unitful.PressureUnits

    coefficients::Array{T3, 3}
    cross_section_unit::Unitful.AreaUnits

end


function show(io::IO, ::MIME"text/plain", absco::ABSCOSpectroscopy4D)

    println(io, "4D ABSCO from file: $(absco.file_name)")
    println(io, "Scale factor: $(absco.scale_factor)")

    println(io, "Pressures: $(length(absco.pressures)) ")
    println(io, "Temperatures: $(size(absco.temperatures, 1))")
    println(io, "Broadener VMRs: $(length(absco.broadener_vmrs))")
    if absco.ww_unit isa Unitful.LengthUnits
        println(io, "Wavelengths: $(length(absco.ww))")
    else
        println(io, "Wavenumbers: $(length(absco.ww))")
    end

end

function show(io::IO, absco::ABSCOSpectroscopy4D)

    println(io, "4D ABSCO from file: $(absco.file_name)")
    println(io, "Scale factor: $(absco.scale_factor)")

    println(io, "Pressures: $(length(absco.pressures)) ")
    println(io, "Temperatures: $(size(absco.temperatures, 1))")
    println(io, "Broadener VMRs: $(length(absco.broadener_vmrs))")
    if absco.ww_unit isa Unitful.LengthUnits
        println(io, "Wavelengths: $(length(absco.ww))")
    else
        println(io, "Wavenumbers: $(length(absco.ww))")
    end

end


function show(io::IO, ::MIME"text/plain", absco::ABSCOSpectroscopy3D)

    println(io, "3D ABSCO from file: $(absco.file_name)")
    println(io, "Scale factor: $(absco.scale_factor)")

    println(io, "Pressures: $(length(absco.pressures)) ")
    println(io, "Temperatures: $(size(absco.temperatures, 1))")
    if absco.ww_unit isa Unitful.LengthUnits
        println(io, "Wavelengths: $(length(absco.ww))")
    else
        println(io, "Wavenumbers: $(length(absco.ww))")
    end

end

function show(io::IO, absco::ABSCOSpectroscopy3D)

    println(io, "3D ABSCO from file: $(absco.file_name)")
    println(io, "Scale factor: $(absco.scale_factor)")

    println(io, "Pressures: $(length(absco.pressures)) ")
    println(io, "Temperatures: $(size(absco.temperatures, 1))")
    if absco.ww_unit isa Unitful.LengthUnits
        println(io, "Wavelengths: $(length(absco.ww))")
    else
        println(io, "Wavenumbers: $(length(absco.ww))")
    end

end