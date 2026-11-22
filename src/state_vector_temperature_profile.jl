"""
$(TYPEDSIGNATURES)

Returns the name of this temperature profile state vector element as a string.
"""
function get_name(sve::TemperatureProfileSVE)
    return "TemperatureProfileSVE ($(sve.level))"
end

"""
$(TYPEDSIGNATURES)

Pretty printing for TemperatureProfileSVE types
"""
function show(io::IO, ::MIME"text/plain", sve::TemperatureProfileSVE)

    println(io, "TemperatureProfileSVE")
    println(io, "Level:             $(sve.level)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    #println(io, "Prior covariance:  $(sve.prior_covariance)")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for TemperatureProfileSVE
"""
function show(io::IO, sve::TemperatureProfileSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end



"""
$(TYPEDSIGNATURES)

Checks whether the atmosphere object `atm` is compatible
with retrieving a temperature profile. The number of levels in both retrieval
and meteorological grid must be the same, and their pressures too!
"""
function check_Tprofile_sanity(atm::AbstractAtmosphere)

    @assert atm.N_level == atm.N_met_level "[TemperatureProfileSVE] Number of levels must be \
    the same for both retrieval grid (N=$(atm.N_level)) and meteorological
    grid (N=$(atm.N_met_level))!"

    return true
    
end

"""
$(TYPEDSIGNATURES)

Returns the positional indices of state vector elements that are a `TemperatureProfileSVE`.
"""
function idx_for_Tprofile_sve(
    sv::AbstractStateVector
)

    # List to be returned
    idx_list = Integer[]
    # Keep a running list of levels, if any appear twice
    # then we have a problem and can warn the user!
    level_list = Integer[]

    for (idx, sve) in StateVectorIterator(sv, TemperatureProfileSVE)
        push!(idx_list, idx)
        push!(level_list, sve.level)
    end

    if length(unique(level_list)) != length(level_list)
        @warn "[INV] There may be duplicate T profile! Following level assignments were found:"
        @warn "[INV]: $(level_list)"
    end

    return idx_list

end

"""
$(TYPEDSIGNATURES)

Returns `N` different `TemperatureProfileSVE` state vector elements, one for each level of
the retrieval atmosphere.
"""
function TemperatureProfileSVE(
    atm::AbstractAtmosphere,
    unit::Unitful.TemperatureUnits;
    type::DataType=Float64
    )

    # This is still experimental: warn the user
    @warn "[SV] Temperature profile retrieval is still in testing phase!"

    # First sanity check happens here
    # (throws error if not fulfilled!)
    check_Tprofile_sanity(atm)

    # Number of levels
    N = atm.N_level

    return [
        TemperatureProfileSVE(
            l,
            unit,
            zero(type),
            zero(type),
            zero(type),
            [zero(type)]
        )
    for l in 1:N]

end


"""
$(TYPEDSIGNATURES)

Returns whether the Jacobian related to `TemperatureProfileSVE` types should be
calculated before convolution happens. Returns `true`.
"""
calculate_jacobian_before_isrf(sve::TemperatureProfileSVE) = true