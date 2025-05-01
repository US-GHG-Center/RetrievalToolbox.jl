function read_hdf5_datasets(
    h5name::String;
    must_contain=nothing,
    must_not_contain=nothing
)


    # Recursively walk through the HDF5 file

    function walk_hdf(obj, d)

        for this_obj in obj

            this_name = split(HDF5.name(this_obj), "/", keepempty=false)[end]

            if this_obj isa HDF5.Group
                d[this_name] = Dict()
                walk_hdf(this_obj, d[this_name])
            end

            if this_obj isa HDF5.Dataset

                if !isnothing(must_contain)
                    # group or dataset must contain "must_contain"
                    if !occursin(must_contain, this_name)
                        continue
                    end
                end

                if !isnothing(must_not_contain)
                    # group or dataset must NOT contain "must_not_contain"
                    if occursin(must_not_contain, this_name)
                        continue
                    end
                end

                d[this_name] = read(this_obj)

            end

        end

    end

    # Recursively walk through dictionary, get rid of empty
    # sub-dicts
    function walk_dict(d)

        for key in keys(d)
            if d[key] isa Dict
                if isempty(d[key])
                    delete!(d, key)
                else
                    walk_dict(d[key])
                end
            end
        end

    end


    h5file = h5open(h5name, "r")
    @debug "Opened HDF5 file: $(h5name)"

    d = Dict()
    walk_hdf(h5file, d)
    close(h5file)

    walk_dict(d)

    return d

end

