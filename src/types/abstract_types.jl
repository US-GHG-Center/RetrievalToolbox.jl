abstract type AbstractSpectroscopy end
abstract type ABSCOSpectroscopy <: AbstractSpectroscopy end

abstract type AbstractAtmosphereElement end
abstract type AbstractAerosolType <: AbstractAtmosphereElement end
abstract type AbstractRayleighScattering <: AbstractAtmosphereElement end

abstract type AbstractBuffer end
abstract type AbstractAtmosphereBuffer <: AbstractBuffer end
abstract type AbstractRTBuffer <: AbstractBuffer end

abstract type AbstractAerosolProperty end

abstract type AbstractSolarModel end

abstract type AbstractSolarSpectrum end

abstract type AbstractSurface end
abstract type BRDFKernel end

abstract type AbstractLocation end

abstract type AbstractAtmosphere end

abstract type AbstractScene end

abstract type AtmosphereScene end

abstract type NoAtmosphereScene end

abstract type AbstractObserver end

abstract type AbstractSpectralWindow end

abstract type AbstractOpticalProperties end

abstract type AbstractISRF end

abstract type AbstractRTMethod end

abstract type AbstractStateVectorElement end

abstract type AbstractStateVector end

abstract type AbstractDispersion end

abstract type AbstractSolver end