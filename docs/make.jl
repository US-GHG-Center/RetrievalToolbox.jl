using Pkg

if isfile(joinpath(@__DIR__, "..", "Project.toml"))
    # Local development - use the local package
    Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
else
    # CI/remote
    # The Github action via the workflow .github/workflows/documentation.yml
    # takes care of adding the module. No need to add anything.
end


using Documenter, RetrievalToolbox
const RE = RetrievalToolbox

ENV["GKSwstype"] = "100" # Inform Plots.jl about headless mode.

DocMeta.setdocmeta!(
    RetrievalToolbox,
    :DocTestSetup,
    :(using RetrievalToolbox);
    recursive=true
)

makedocs(
    sitename="RetrievalToolbox.jl",
    remotes=nothing,
    #modules=[RetrievalToolbox],
    doctest=true,
    pages = [
        "Main" => "index.md",
        "Concepts" =>
            [
                "Fundamentals" => joinpath("concepts", "fundamentals.md"),
                "Core Concepts" => joinpath("concepts", "core_concepts.md"),
                "Radiance" => joinpath("concepts", "radiance.md"),
                "Scattering Phasefunction" => joinpath("concepts" , "phasefunction.md"),
            ],
        "Design" => joinpath("design", "design.md"),
        "Functions" =>
            [
                "State Vector Functions" => joinpath("functions", "state_vector_functions.md"),
                "Atmosphere Functions" => joinpath("functions", "atmosphere_functions.md"),
                "Instrument Functions" => joinpath("functions", "instrument_functions.md")
            ],
        "Types" =>
            [
                "Atmosphere Types" => joinpath("types", "atmosphere_types.md"),
                "Buffer Types" => joinpath("types", "buffer_types.md"),
                "Dispersion Types" => joinpath("types", "dispersion_types.md"),
                "Radiative Transfer Method Types" => joinpath("types", "RT_types.md"),
                "State Vector Types" => joinpath("types", "state_vector_types.md"),
                "Surface Types" => joinpath("types", "surface_types.md"),
            ],
        "Working with Julia" =>
            [
                "Develop" => joinpath("julia", "develop.md"),
                "Dictionaries" => joinpath("julia", "dicts.md"),
                "Units" => joinpath("julia", "units.md"),
            ],
        "Pitfalls" =>
            [
                "Pitfalls" => joinpath("pitfalls", "pitfalls.md"),
            ],
    ]
)

deploydocs(;
    repo = "github.com/US-GHG-Center/RetrievalToolbox.jl.git",
    versions = ["stable" => "v^", "v#.#", "main" => "main", "dev" => "dev"],
    push_preview = false,
)

# Remove package from Project.toml
Pkg.rm("RetrievalToolbox")
