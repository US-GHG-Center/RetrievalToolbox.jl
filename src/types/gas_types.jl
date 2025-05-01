"""
Gas absorber type. Requires an associated spectroscopy object
along with pressure levels and the corresponding volume mixing ratios
on each of those levels.

$(TYPEDFIELDS)
"""
struct GasAbsorber{T} <: AbstractAtmosphereElement where T

    "Name of the gas"
    gas_name::String
    "Spectroscopy object that contains the cross sections for this gas"
    spectroscopy::AbstractSpectroscopy
    "VMR level profile"
    vmr_levels::Vector{T}
    "VMR units"
    vmr_unit::Unitful.DimensionlessUnits

end



"""
Pretty printing for gas absorber types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", gas::GasAbsorber)

    println(io, "Gas absorber: $(gas.gas_name)")
    if hasproperty(gas.spectroscopy, :file_name)
        println(io, "Associated spectroscopy: $(gas.spectroscopy.file_name)")
    else
        println(io, "Associated spectroscopy: $(gas.spectroscopy)")
    end

end

"""
Brief pretty printing for gas absorber types

$(SIGNATURES)
"""
function show(io::IO, gas::GasAbsorber)

    print(io, "GasAbsorber: $(gas.gas_name)")

end
