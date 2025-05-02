"""
$(TYPEDSIGNATURES)

Generic function to apply a radiance correction to an at-instrument spectrum, for example
stray light or a zero-level offset. The generic implementation does nothing, and users
should implement their own in accordance with the state vector element `sve`.
"""
function apply_radiance_correction!(
    rt_buf::AbstractRTBuffer,
    sve::AbstractStateVectorElement
    )

    # Default implementation, which does nothing
    return nothing

end



"""
$(TYPEDSIGNATURES)

Radiance "correction" due to a zero-level offset as given by the
`ZeroLevelOffsetPolynomialSVE` object. The spectrally-dependent radiance is added to the
total radiance, and the intensity component only.
"""
function apply_radiance_correction!(
    rt_buf::AbstractRTBuffer,
    sve::ZeroLevelOffsetPolynomialSVE
    )

    # Perform the ZLO calculation
    calculate_zlo!(
        rt_buf.radiance.I,
        rt_buf.radiance_unit,
        sve,
        rt_buf.dispersion[sve.swin]
    )

end