# Instrument functions

RetrievalToolbox is mostly instrument-agnostic, meaning that none of the functions are
written for specific instruments, even though a certain type of instrument class is
assumed for the time being.

Nevertheless, in order to create a fully functioning retrieval algorithm, the forward
model must be capable of producing instrument-level radiances. The functions listed below
provide ways of doing so, assuming the high-resolution radiances and Jacobians have been
calculated already.

Note that some of the functions below may refer to quantities whose meaning should be
understood in detail, for further reading consult [`Core Algorithm Concepts`](@ref),
especially the section on *Pixels, Spectral Samples and Dispersion*.

## List of instrument functions
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["instrument.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```