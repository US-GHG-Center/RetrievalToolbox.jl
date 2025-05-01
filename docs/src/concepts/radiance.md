# Working with radiance types

The way RetrievalToolbox handles radiances is different, when compared to other objects such as


There are, at the moment, two radiance types that users can utilize to store calculated or measured radiance:

```julia
struct ScalarRadiance{T<:AbstractFloat} <: Radiance

    I::Vector{T}

end
```

and

```julia
struct VectorRadiance{T<:AbstractFloat} <: Radiance

    I::Vector{T}
    Q::Vector{T}
    U::Vector{T}

end
```

Depending on the situation, users can work with `ScalarRadiance` when polarization does not need to be accounted for, or use the `VectorRadiance` type if polarization is needed. Note that as opposed to many other types, these radiance containers do not have a unit field, but that might change in future versions.

The main idea is to use whichever type is most natural for the specific application. While it would have been possible to only use a single radiance type which considers all components of the Stokes vector, it would require users to pay attention to Stokes components *Q* and *U* without ever truly needing them.



## Basic examples

For a very basic example, let us first generate some scalar radiance object using the internal constructor (not included in the above type definitions).

```@repl rad
using RetrievalToolbox # hide
const RE = RetrievalToolbox # hide
T = Float64;
N = 4;
s = RE.ScalarRadiance(T, N)
```

`s` is now our radiance object, which only has one field, namely `I`, representing the intensity for some number of spectral indices.

Similarly, we can create a vector radiance that represents the first three components of the Stokes vector:

```@repl rad; continued = true
v = RE.VectorRadiance(T, N)
```

Functions that manipulate radiance can easily access the Stokes components inside either `s` or `v` via the explicit dot syntax.

```@repl rad; continued = true
s.I[:] .= 0.1;
v.I[:] .= 0.2;
v.Q[:] .= -0.3;
v.U[:] .= 0.03;
```



As with many types of RetrievalToolbox, the underlying fields are simple vectors, such that highly performant operations can be done with them. For now, RetrievalToolbox makes use of `LoopVectorization.jl`, which speeds up certain looped computations dramatically, and we can easily use e.g. the `@turbo` macro to wrap such loops. Below example shows a simple in-place calculation involving both `s` and `v` in a single loop body.

```@repl rad; continued = true
using LoopVectorization

@turbo for i in eachindex(s.I)
    s.I[i] += v.I[i]
    v.U[i] /= 2
    v.Q[i] /= 2
end
s
v
```


## Using `StokesIterator`