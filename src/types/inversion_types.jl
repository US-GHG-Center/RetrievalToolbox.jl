include("inversion_IMAP_types.jl")
include("inversion_LM_types.jl")

"""
A type to hold matrices and other vectors for
Rodgers-type OE error analysis.

$(TYPEDFIELDS)

Note that these fields are kept deliberately as
Abstract-type quantities, such that views on arrays
and vectors can be used.

"""
struct OEQuantities

    K::AbstractMatrix # Jacobian
    Se::AbstractMatrix # Error covariance
    Sa::AbstractMatrix # Prior covariance
    Shat::AbstractMatrix # Posterior covariance
    G::AbstractMatrix # Gain matrix
    AK::AbstractMatrix # Averaging kernel matrix
    SV_ucert::AbstractVector # Posterior SV uncertainty
    SV_ucert_noise::AbstractVector # SV uncertainty due to noise
    SV_ucert_smoothing::AbstractVector # SV uncertainty due to smoothing

end
