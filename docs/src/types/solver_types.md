# [Solver Types](@id solver_types)

To facilitate inversions of measured radiances in an organized way, RetrievalToolbox provides dedicated types that implement solver methods. These solvers can be understood as light-weight wrapper methods that merely implement the needed linear algebra to obtain optimized state vectors, and are not strongly tied to a specific forward model function, or even a specific way how that forward model function may be written.

Users may write their own solver types to implement new inversion methods. In the type hierarchy tree, they should be a subtype of `AbstractSolver`. In order to allow for seamless functionality with some common functions from [Inversion Functions](@ref), we suggest to model the new solver type according to the `IMAPSolver` below, however advanced users may be able move away from that template. Once the type definition is written, users should then write a number of functions that implement the calculation of state vector updates, error analysis and so on. We recommend looking at the `src/inversion_IMAP.jl` file to get a feeling of how these functions should be written.

Forward model functions can be rather complex, we highly recommend to study existing example implementations available on the [RetrievalToolbox Org Github](https://github.com/RetrievalToolbox).

## The forward model function

Most solver types should have a function `forward_model` as a field within the solver type. These functions are supposed to have only one single argument: the state vector (of type `<:AbstractStateVector`). The function `forward_model` may have any number of **additional keyword arguments**, which in Julia are distinguished by a semi-colon. Further, it is expected by some other functions that the forward model return `true` if it executed successfully. The reasoning behind this is that forward model functions may be highly complex and it is very helpful to signal to the inversion methods that the forward model failed.

For example, a forward model function may look like this:

```julia
function forward_model!(
    sv::RetrievalStatevector;
    arg1::EarthAtmosphereBuffer,
    arg2::Bool,
    arg3::Float64
)

    # user code ..

    return true
end
```

Note the exclamation point in the function name, which signifies (by Julia convention) that the function modifies at least one of its arguments. When creating a solver object, one would then pass only the function name:

```julia
# .. set up

solver = IMAPSolver(
    forward_model!,
    sv,
    .. # other arguments ignored here for brevity
)
```

## Iterating and handling keyword arguments

Solver types (inversion methods) that are part of RetrievalToolbox come with a `next_iteration!` function that is implemented for each solver type individually. Those are really the core functions that call the forward model, and then perform the calculations to update the state vector elements.

We generally suggest following way of organizing keyword arguments to forward functions as well as perform an iteration. First, establish the keyword arguments through a named tuple which uses exactly the same names as the arguments in the forward model function:

```julia
fm_kwargs = (
    arg1 = my_buffer,
    arg2 = True,
    arg3 = 10.5
)
```

Then call the routine to perform the next step of the inversion:

```julia
next_iteration!(solver; fm_kwargs)
```

**Important!** The `next_iteration!` function does **not** take the specific keyword arguments like the forward model function does, instead it takes a (named) tuple. Internally, `next_iteration!` then expands the named tuple to call the forward model function with the specific keyword arguments. If the forward model needs to be called for another reason, outside of an iteration step, the keyword arguments must be expanded via the so-called splat operator `...`:

```julia
forward_model!(sv; fm_kwargs...)
```

or (equivalently)

```julia
solver.forward_model(sv; fm_kwargs...)
```



## Iterative maximum a-posteriori (IMAP) solver

```@docs
RE.IMAPSolver
```

```@autodocs
Modules = [RetrievalToolbox]
Pages = ["inversion_IMAP.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```

## Levenberg-Marquardt solver

```@docs
RE.LMSolver
```

```@autodocs
Modules = [RetrievalToolbox]
Pages = ["inversion_LM.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```
