struct RayleighScattering <: AbstractRayleighScattering

    # No further parameters needed

end


"""
Pretty printing for RayleighScattering types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", ray::RayleighScattering)

    println(io, "RayleighScattering")

end

"""
Brief pretty printing for gas absorber types

$(SIGNATURES)
"""
function show(io::IO, ray::RayleighScattering)

    print(io, "RayleighScattering")

end