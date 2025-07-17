using Documenter, RetrievalToolbox
const RE = RetrievalToolbox


makedocs(
    sitename="RetrievalToolbox.jl",
    remotes=nothing,
    pages = [
        "Main" => "index.md",
        "Design" => joinpath("design", "design.md"),
        "Concepts" =>
            [
                "Core Concepts" => joinpath("concepts", "core_concepts.md"),
                "Radiance" => joinpath("concepts", "radiance.md"),
                "Scattering Phasefunction" => joinpath("concepts" , "phasefunction.md"),
            ],
        "Functions" =>
            [
                "State Vector Functions" => joinpath("functions", "state_vector_functions.md"),
            ],
        "Types" =>
            [
                "Buffer Types" => joinpath("types", "buffer_types.md"),
                "State Vector Types" => joinpath("types", "state_vector_types.md"),
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
    versions = ["stable" => "v^", "v#.#", "dev" => "main"],
    push_preview = false,
)