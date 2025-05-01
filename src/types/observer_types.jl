"""
A type to hold a satellite observer, including its
viewing angles, position and velocity vectors.

$(TYPEDEF)

$(TYPEDFIELDS)

"""
mutable struct SatelliteObserver{T1, T2} <: AbstractObserver

    viewing_zenith::T1
    viewing_azimuth::T1
    satellite_velocity::Vector{T2}
    satellite_position::Vector{T2}

end

struct UplookingGroundObserver <: AbstractObserver end

"""
A type to hold an observer that sits somewhere at a
finite height in the atmosphere (e.g. an aircraft).
This type holds the viewing angles, as well as the
observer pressure and altitude, plus the velocity vector.

$(TYPEDEF)

$(TYPEDFIELDS)

"""
mutable struct AtHeightObserver{T} <: AbstractObserver where T <: AbstractFloat

    viewing_zenith::T
    viewing_azimuth::T
    observer_altitude::T
    observer_pressure::T
    observer_velocity::Vector{T}

end
