#=
    See https://doi.org/10.1364/AO.44.003320
    or https://opg.optica.org/ao/viewmedia.cfm?uri=ao-44-16-3320
    (Tomasi et al. 2005)
    for many of these calculations.
=#

"""
$(TYPEDSIGNATURES)

Computes the King factor according to Tomasi et al. 2005,
`https://doi.org/10.1364/AO.44.003320`. Requires `wl_input` as a length-valued input.
"""
function King_factor(wl_input::Unitful.Length)

    # This ignores water vapor

    # Convert into microns to be consistent with the equations in Tomasi et al.
    wl = ustrip(u"µm", wl_input)

    N2_C = 78.084
    O2_C = 20.946
    Ar_C = 1.0
    CO2_C = 0.0004

    F_Ar = 0.934
    F_CO2 = 1.15

    wl2 = @. wl * wl # squared wavelength [um^2]
    wl4 = @. wl2 * wl2 # wavelength^4 [um^4]

    F_N2 = @. 1.034 + 3.17 * 0.0001 / wl2
    F_O2 = @. 1.096 + 1.385 * 0.001 / wl2 + 1.448 * 0.0001 / wl4

    fac = @. N2_C * F_N2 + O2_C * F_O2 + Ar_C * F_Ar + CO2_C * F_CO2
    fac = @. fac / (N2_C + O2_C + Ar_C + CO2_C)

    # The returned quantity is of unit 1.
    return fac

end


function calculate_rayleigh_depolf(wl::Unitful.Length)

    # Wavelength here can be any unit, since the King_factor function takes care of the
    # conversion.
    F = King_factor(wl)
    depol_r = (6.0 * F - 6.0) / (7.0 * F + 3.0)
    return depol_r

end


function calculate_rayleigh_sigma(wn::Unitful.Wavenumber)

    # Convert to wavelength
    wn = 1 / wn |> u"µm"
    # Calculate in wavelength space
    return calculate_rayleigh_sigma(wn)

end



function calculate_rayleigh_sigma(
    wl::Unitful.Length # Wavelength
    )

    Nair = 2.546899e19 #u"1/cm^3"

    pi3 = pi * pi * pi

    F = King_factor(wl)

    # Wavelength must be in units of micro meters
    # for the calculation of the refractive index

    wl2 = ustrip(u"µm", wl)^2
    wli2 = 1 / wl2

    # Refractive index (with unit 1)
    # This is equation (9) in Tomasi et al., and "wli2",
    # is 1\lambda^2 where lambda must be in micro meters.
    n300_1 = 8_060.51 + 2_480_990.0 / (132.274 - wli2) +
        17_455.7 / (39.32957 - wli2)
    n300_1 = n300_1 * 1.0e-8
    ns = 1.0 + n300_1

    # In order for the result to be in cm^2/molecule,
    # we must convert the wavelength into cm.
    wl4_cm = ustrip(u"cm", wl)^4

    ρ = calculate_rayleigh_depolf(wl)
    ray_sigma = (24 * pi3 * (ns*ns - 1) * (ns*ns - 1)) /
        (wl4_cm * Nair * Nair * (ns*ns + 2) * (ns*ns + 2)) * (6 + 3 * ρ) / (6 - 7 * ρ)

    # The returned quantity is in cm^2 per molecule
    # (drop per molecule)
    return ray_sigma*u"cm^2"

end


#=
    NOTE
    ====

    Below functions are "doubling up" a little and define Rayleigh calculations for both
    wavelength and wavenumbers. This is mostly for performance reasons; it guarantees
    type stability for the wavenumber of wavelength quantity. One could write one function
    only that checks the type inside the function, but this would make the function
    likely slower. So for the sake of performance, we have wrapper functions for both
    wavelength and wavenumber spectral units. HOWEVER there is only *one* function that
    actually implements the Rayleigh coefficient calculations, so there is one place
    where mistakes can occur there..

=#




function create_rayleigh_coefs(
    wn::Unitful.Wavenumber,
    polarized::Bool
)

    if polarized
        coef = zeros(3, 6)
    else
        coef = zeros(3, 1)
    end

    create_rayleigh_coefs!(coef, wn, polarized)

    return coef

end


function create_rayleigh_coefs(
    wl::Unitful.Length,
    polarized::Bool
)

    if polarized
        coef = zeros(3, 6)
    else
        coef = zeros(3, 1)
    end

    create_rayleigh_coefs!(coef, wl, polarized)

    return coef

end


function create_rayleigh_coefs!(
    coef::Matrix,
    wn::Unitful.Wavenumber,
    polarized::Bool
    )

    wl = 1 / wn |> u"µm"
    create_rayleigh_coefs!(coef, wl, polarized)

end


"""
$(TYPEDSIGNATURES)

Calculates the Rayleigh scattering matrix into an existing array of size (3,1) or (3,6),
depending on whether the user needs polarization. This modifies the matrix `coef`
in-place.
"""
function create_rayleigh_coefs!(
    coef::Matrix,
    wl::Unitful.Length,
    polarized::Bool,
    )

    # Zero out!
    @views coef[:,:] .= 0
    # Get the depolarization factor for the requested wavelength
    depolf = calculate_rayleigh_depolf(wl)

    #= Indices: [moment, (Greek) element]
        Greek elements must be ordered:
        β (l,1), α (l,2), ζ (l,3), δ (l,4), -γ (l,5), -ε (l,6)

        Be mindful of the fact that we use 1-based indexing here, whereas most literature
        will likely use 0-based indexing for the moments `l`, thus β0 in some paper will
        probably correspond to coef[1,1], or δ8 correspond to coef[9,4].

        As further explained in the documentation, we also observe the Siewert convention
        regarding the construction of the phase matrix expansion (for the time being, due
        to compatibility with XRTM). The elements (l,5) and (l,6) **must** correspond to
        -γ_l and -ε_l, and include the minus sign!

        Note the source from RtRetrievalFramework/Support/rayleigh_greek_moment.cc, and
        the VLIDORT documentation section 3.1.1.
    =#

    coef[1,1] = 1.0 # β0
    coef[2,1] = 1e-11 # β1  Apparently this might be needed to keep RT from breaking..
    coef[3,1] = (1.0 - depolf) / (2.0 + depolf) # β2

    if polarized
        coef[3,2] = 6.0 * (1.0 - depolf) / (2.0 + depolf) # α2
        coef[2,4] = 3.0 * (1.0 - 2.0 * depolf) / (2.0 + depolf) # δ1
        coef[3,5] = sqrt(6.0) * (1.0 - depolf) / (2.0 + depolf) # -γ2
        # Note the lack of minus in front!
        # XRTM uses requires "-γ" as input, rather than "γ"
    end

end