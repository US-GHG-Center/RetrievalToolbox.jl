#=
This script produces the waveform file from the data published in Magney et al. 2019,
https://doi.org/10.1029/2019JG005029

The underlying files are the same as found on the data repository at
https://doi.org/10.22002/D1.1226

Users of RetrievalToolbox should have no need to run this script, since the output is
already incorporated statically into the code. Those who still want to re-produce the
waveform file must first download the four "FluorescenceSpectra*.csv" files from the
repository above and then run this code through the `do_all()` function.
=#

using CSV
using DataFrames
using DocStringExtensions
using Glob
using LinearAlgebra
using Statistics

"""
$(TYPEDSIGNATURES)

Creates a `DataFrame` from the Magney et al. data file
"""
function create_dataframe(fname::String)

    # Grab data and join into one long string
    raw = readlines(fname)
    # Get rid of \ufeff characters
    raw_filtered = replace.(raw, Ref('\ufeff' => ""))
    # Drop samples that have just a "\"
    bad_rows = Int[]
    for (i,line) in enumerate(raw_filtered)
        if occursin("\\", line)
            push!(bad_rows, i)
        end
    end

    # Get rid of bad rows
    deleteat!(raw_filtered, bad_rows)

    return CSV.File(IOBuffer(join(raw_filtered, '\n'))) |> DataFrame
end

"""
$(TYPEDSIGNATURES)

Creates the matrix of samples vs. radiance from a `DataFrame`. After transposition,
columns will indicate samples, and rows will indicate datapoints (radiance at wavelength).
"""
function create_matrix(df)

    # Find out which wavelength columns have successfully been
    # parsed into Floats.
    # Note that we transpose the matrix here so we can go on to
    # calculate the PCA.
    wl_c = names(df, Float64)
    return df[2:end, wl_c] |> Array |> transpose

end

"""
$(TYPEDSIGNATURES)

Normalizes the matrix by mean-removal and division by standard deviation along columns.
Returns the new normalized matrix.
"""
function normalize_matrix(m)
    return (m .- mean(m, dims=1)) ./ std(m, dims=1)
end

"""
$(TYPEDSIGNATURES)

Gets the vector of wavelengths corresponding to the samples for a `DataFrame`.
"""
function get_wavelengths(df)
    # Find out which wavelength columns have successfully been
    # parsed into Floats
    wl_c = names(df, Float64)
    return df[1, wl_c] |> Vector
end

"""
$(TYPEDSIGNATURES)

Perform the PCA on a normalized matrix
"""
function create_waveform_through_PCA(mbar)

    eig = eigen(cov(mbar), sortby=-)
    proj = mbar * eig.vectors

    @info "PC1 variance explained: $(eig.values[1] / sum(eig.values))"

    # Grab first PC
    PC1 = proj[:,1]

    # Flip, Offset and Normalize
    if PC1[argmax(abs.(PC1))] < 0
        PC1[:] .*= -1
    end

    PC1[:] .+= abs.(minimum(PC1))
    PC1[:] ./= maximum(PC1)

    return PC1

end

"""
Perform the full operation to produce the "SIF_waveform.csv" file.
"""
function do_all()

    fnames = Glob.glob("Fluorescence*.csv")

    dfs = fnames .|> create_dataframe
    wls = dfs .|> get_wavelengths

    # Make sure all wavelengths are the same
    if !(all(diff(hcat(wls...), dims=2) .== 0))
        return nothing
    end

    ms = dfs .|> create_matrix
    msbar = hcat(ms...) |> normalize_matrix

    PC1 = msbar |> create_waveform_through_PCA

    open("SIF_waveform.csv", "w") do fout
        for i in 1:length(wl)
            println(fout, "$(wl[i]),$(PC1[i])")
        end
    end

end