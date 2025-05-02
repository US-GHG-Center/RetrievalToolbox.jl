"""
Returns the name of this dispersion state vector element as a string.

$(SIGNATURES)

"""
function get_name(sve::GasVMRProfileSVE)
    return "GasVMRProfileSVE $(sve.gas.gas_name) ($(sve.level))"
end

"""
Pretty printing for GasVMRProfileSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::GasVMRProfileSVE)

    println(io, "GasVMRProfileSVE")
    println(io, "Level:             $(sve.level)")
    println(io, "Corresponding gas: $(sve.gas.gas_name)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    #println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
Brief pretty printing for GasVMRProfileSVE

$(SIGNATURES)
"""
function show(io::IO, sve::GasVMRProfileSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

"""
Returns `N` different `GasVMRProfileSVE` state vector elements, one for each level of
the retrieval atmosphere.
"""
function GasVMRProfileSVE(
    N::Integer,
    gas::GasAbsorber,
    unit::Unitful.DimensionlessUnits;
    type::DataType=Float64
    )

    return [
        GasVMRProfileSVE(
            l,
            gas,
            unit,
            zero(type),
            zero(type),
            zero(type),
            [zero(type)]
        )
    for l in 1:N]

end

"""

"""
function calculate_jacobian_before_isrf(sve::GasVMRProfileSVE)

    return true

end


"""
Returns the positional indices of state vector elements that are a `GasVMRProfileSVE` and
belong to a `GasAbsorber` `gas`. This allows easy retrieval of the VMR profile SVEs for
some gas
"""
function idx_for_profile_sve(
    gas::GasAbsorber,
    sv::AbstractStateVector
)

    idx_list = Integer[]
    for (idx,sve) in StateVectorIterator(sv, GasVMRProfileSVE)

        if sve.gas === gas
            push!(idx_list, idx)
        end

    end

    return idx_list

end