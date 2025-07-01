#=
    This section of the module dynamically generates new Base.getproperty functions, such
    that users can access appropriate struct fields as .wavelength or .wavenumber,
    depending on the unit chosen during instantiation.

    This allows a "natural" treatment of spectral quantities in both wavelength and
    waveumber spaces without having to resort to some placeholder name.

    IMPORTANT Do not use s.x or getproperty(s, x) within the new, overloaded function,
    since it will cause an infinite recursive call to getproperty.
=#

# Loop through every struct within the module
for name in names(@__MODULE__, all=true)

    if getfield(@__MODULE__, name) isa Type

        # Abstract types cannot be processed this way so, skip those
        try (fieldcount(getfield(@__MODULE__, name)))
            @debug "[MISC] $(name) will be checked for magical ww."
        catch
            @debug "[MISC] $(name) is not compatible here."
            continue
        end

        # Check if this type has a "ww" in it
        has_ww = false
        for fname_sym in fieldnames(getfield(@__MODULE__, name))
            if startswith(String(fname_sym), "ww")
                has_ww = true
                break
            end
        end

        # Skip the rest if no "ww" is found in any of the field names
        if !has_ww
            continue
        end

        # For each type that has a "ww" field in it, we need
        # to write one single Base.getproperty function.

        function_string = """
            function Base.getproperty(s::$(name), x::Symbol)

        """

        for fname_sym in fieldnames(getfield(@__MODULE__, name))
            fname_str = String(fname_sym)

            if startswith(fname_str, "ww")

                nl1 = replace(fname_str, "ww" => "λ")
                nl2 = replace(fname_str, "ww" => "wavelength")

                nn1 = replace(fname_str, "ww" => "ν")
                nn2 = replace(fname_str, "ww" => "wavenumber")

                @debug "[MISC] Creating λ/ν magic accessors for $(name)::$(fname_str)"
                @debug "[MISC] $(name)::$(fname_str) -> $(name)::$(nl1), $(name)::$(nn1)"
                @debug "[MISC] $(name)::$(fname_str) -> $(name)::$(nl2), $(name)::$(nn2)"

                # short-circuit evaluation
                function_string *= """
                        ((x === :$(nl1)) || (x === :$(nl2))) &&
                            getfield(s, :ww_unit) isa Unitful.LengthUnits &&
                            return getfield(s, :$(fname_str))
                        ((x === :$(nn1)) || (x === :$(nn2))) &&
                            getfield(s, :ww_unit) isa Unitful.WavenumberUnits &&
                            return getfield(s, :$(fname_str))
                        ((x === :$(nl1)) || (x === :$(nl2))) &&
                            getfield(s, :ww_unit) isa Unitful.WavenumberUnits &&
                            error("Wrong unit! Expected wavenumber, but got λ.")
                        ((x === :$(nn1)) || (x === :$(nn2))) &&
                            getfield(s, :ww_unit) isa Unitful.LengthUnits &&
                            error("Wrong unit! Expected wavelength, but got ν.")
                """

            end

        end

        function_string *= """

                # default dispatch
                return getfield(s, x)

            end
        """

        # Evaluate the new, overloaded Base.getproperty function!
        eval(Meta.parse(function_string))

    end

end