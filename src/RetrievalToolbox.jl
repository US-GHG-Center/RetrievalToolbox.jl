module RetrievalToolbox

import Base: length, show, iterate, convert, getindex

using AstroLib
using CSV
using Dates
using DocStringExtensions
using HDF5
using Interpolations
using Lazy
using LinearAlgebra
using LoopVectorization
using NCDatasets
using Polynomials
using PrettyTables
using Printf
using ProgressMeter
using QuadGK
using SharedArrays
using Statistics
using StaticArrays
using Unitful


#=
    Make the "photons" unit available outside of
    the module scope, otherwise users will have to
    use RetrievalToolbox.ph instead of u"ph".
=#
export ph

# Define spectral radiance units
@unit ph "ph" Photons 1u"1" true;
Unitful.register(RetrievalToolbox)

# Spectral radiance
rad_W = u"W/m^2/sr/µm"
rad_ph = u"ph/s/m^2/sr/µm"

# Spectral irradiance
irrad_W = u"W/m^2/µm"
irrad_ph = u"ph/s/m^2/µm"

# Definite unit types analogous to e.g. Unitful.LengthUnits

Unitful.register(@__MODULE__)

# Function to filter out methods that are of no interest in the documentatiton
function DocFilter(f)

    fname = String(nameof(f))

    # No need to list `show` functions for printing
    fname == "show" && return false
    # No need to show internals (_do_this_and_that)
    startswith(fname, "_") && return false

    # Otherwise simply return true!
    return true

end


function __init__()

    # Register photon units at runtime
    Unitful.register(@__MODULE__)

    #=
        If the XRTM radiative transfer library is available, we can load it here. The
        location of the XRTM interface is determined via the environment variable
        XRTM_PATH. Users must set that variable to either the location of the XRTM root
        path, or directly the path which contains the interface XRTM.jl *and* the
        compiled and linked library.
    =#

    have_XRTM = false

    if !haskey(ENV, "XRTM_PATH")

        @debug "XRTM_PATH environment variable was not found! Not loading XRTM!"

    else

        # Try some possibilities where XRTM could be hiding
        xrtm_path_options = String[]

        push!(xrtm_path_options, ENV["XRTM_PATH"])
        push!(xrtm_path_options, joinpath(ENV["XRTM_PATH"], "interfaces"))

        # Loop through options and try to include the XRTM.jl file
        for s in xrtm_path_options

            xrtm_path = joinpath(s, "XRTM.jl")

            if isfile(xrtm_path)
                @debug "Loading XRTM from $(s)"
                include(xrtm_path)

                # Also set the environment variable for this Julia session to find the
                # compiled XRTM library itself.
                push!(Base.DL_LOAD_PATH, s)
                # Must set this global, because it is a module-wide variable unfortunately
                have_XRTM = true
                break
            else
                @debug "[XRTM] No XRTM.jl in $(s)"
            end
        end

        return nothing

    end

    if !(have_XRTM)
        @debug "Could not find XRTM module!"
        @debug "Calls to XRTM library functions will crash the session!"
    end

end




# These will be used by the codes below
include("const.jl") # Physical constants
include("misc.jl") # Miscellaneous routines

### MAIN BODY OF RetrievalToolbox ###

# Definition of abstract types, which
# have to be defined before concrete types.
include("types/abstract_types.jl")

# Definition of concrete types
include("types/gas_types.jl")
include("types/state_vector_types.jl")
include("types/spectral_window_types.jl")
include("types/dispersion_types.jl")
include("types/spectroscopy_types.jl")
include("types/aerosol_types.jl")
include("types/rayleigh_types.jl")
include("types/solar_model_types.jl")
include("types/observer_types.jl")
include("types/surface_types.jl")
include("types/location_types.jl")
include("types/optical_properties_types.jl")
include("types/atmosphere_types.jl")
include("types/scene_types.jl")
include("types/instrument_types.jl")
include("types/radiance_types.jl")
include("types/RT_types.jl")
include("types/sif_types.jl")

include("types/inversion_types.jl")

include("types/buffer_types.jl")


#=
    Here we insert the magic, which allows us to refer to spectral
    domain fields of types with either .wavelength or .wavenumber,
    depending on the unit that was used to instantiate it.
=#
include("magic_ww.jl")

# Implementations of functions

# State vector elements
include("state_vector_surface_albedo.jl")
include("state_vector_dispersion.jl")
include("state_vector_ils_stretch.jl")
include("state_vector_surface_pressure.jl")
include("state_vector_zero_level_offset.jl")
include("state_vector_gas_scale.jl")
include("state_vector_gas_profile.jl")
include("state_vector_solar_scale.jl")
include("state_vector_temperature_offset.jl")
include("state_vector_aerosol_od.jl")
include("state_vector_aerosol_height.jl")
include("state_vector_aerosol_width.jl")
include("state_vector_brdf_polynomial.jl")
include("state_vector_sif_radiance.jl")

include("state_vector.jl")

# Rest
include("spectral_window.jl")
include("dispersion.jl")
include("spectroscopy.jl")
include("gas.jl")
include("aerosol.jl")
include("rayleigh.jl")
include("solar_model.jl")
include("observer.jl")
include("surface.jl")
include("location.jl")
include("scene.jl")
include("optical_properties.jl")
include("atmosphere.jl")
include("instrument.jl")
include("radiance.jl")
include("RT.jl")
include("radiance_correction.jl")
include("SIFRadiance.jl")

include("inversion.jl")

include("buffer.jl")

# Misc other tools that help you run the code
include("hdf5_tools.jl")
include("OCO_tools.jl")


# We export ALL symbols in the module that do not start with an underscore
# (This bit of code is taken from JuMP https://github.com/jump-dev/JuMP.jl)
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end



@info "Loaded RetrievalToolbox.jl, $(Dates.now())"

end # module
