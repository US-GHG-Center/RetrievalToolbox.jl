# Working with radiance types

There are, at the moment, two radiance types that users can utilize to store calculated or measured radiance: `ScalarRadiance` and `VectorRadiance`. Both very similar to generic arrays, in fact they inherit many of the properties and functions from Julia's arrays.

Depending on the situation, users can work with `ScalarRadiance` when polarization does not need to be accounted for, or use the `VectorRadiance` type if polarization is needed. Note that as opposed to many other types, these radiance containers do not have a unit field, but that might change in future versions.

The main idea is to use whichever type is most natural for the specific application. While it would have been possible to only use a single radiance type which considers all components of the Stokes vector, it would require users to pay attention to Stokes components *Q* and *U* without ever truly needing them.

Both `ScalarRadiance` and `VectorRadiance` have only one field, named `S` which itself is the underlying array of some type `T` to store the radiance. For `ScalarRadiance`, `S` is a vector (intensity *I* only), and `VectorRadiance`, `S` is a 3-column array representing the *I*, *Q* and *U* components. Both radiance types possess appropriate accessor functions so users can use the more "natural" I, Q, U notation to access the Stokes components in addition to simply using the type field `.S`.

## Basic examples

For a very basic example, let us first generate some scalar radiance object using the internal constructor.

```@repl rad
using RetrievalToolbox # hide
const RE = RetrievalToolbox # hide
T = Float64;
N = 4;
s = RE.ScalarRadiance(T, N)
```

`s` is now our radiance object, which only has one field, namely `S`, representing the intensity for some number of spectral indices.

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


The main reason why the underlying object `.S` is an Array/Vector-type for both is to facilitate highly performant operations. For now, RetrievalToolbox makes use of `LoopVectorization.jl`, which speeds up certain looped computations dramatically, and we can easily use e.g. the `@turbo` macro to wrap such loops. Having both radiance types look similar under the hood allows users to write fast functions that are agnostic to the type of radiance, and thus work for either. Yet, it is still possible to dispatch on a particular radiance type, so that specific calculations can be done according to whether some radiance object is a `ScalarRadiance` or a `VectorRadiance`

For example, let's write a function that simply doubles all components of a radiance object:

```@repl rad; continued = true
using LoopVectorization

function double_rad!(r::Radiance)
    @turbo for i in axes(r, 1) # Loop over spectral points
        for j in axes(r, 2) # Loop over components (if applicable)
            r.S[i,j] *= 2
        end
    end
end

double_rad!(s)
show(s)
double_rad!(v)
show(v)
```

Note the following: first, we make sure that this function is only to be used on radiance types, so we use the `Radiance` type (which is simply a union between `ScalarRadiance` and `VectorRadiance`) to restrict the function argument `r`. Then, we write a @turbo-accelerated loop that stretches over two dimensions of r. Even for the one-dimensional scalar radiance `s`, a second dimension is accessible as long as that second dimension is accessed by the index `[1]`, trying to access another index will lead to a `BoundsError`. Lastly we do not simply write `r[i,j] *= 2`, but explicitly state the underlying array object `r.S` and write `r.S[i,j] *= 2`. While both versions would work in this case, the `@turbo` macro would not be able to produce accelerated code if we omit the `.S`.

!!! note
    When writing code to perform operations on or with some radiance object `r`, best performance is achieved when explicitly using the underlying array object `r.S`.


The ability to write a two-dimensional loop body even for the one-dimensional scalar radiance is crucial and allows users to write radiance type-agnostic code. Mixing two radiance objects of different radiance type is also easy to do, but requires some more checks on the possible types. Note that **as a conscious design decision, there are currently no arithmetic operations defined on radiance objects of different types**, so users will have to explicitly write those operations. While one can perform e.g. additions on two radiance objects of **the same type**, operations on different types will generally fail.

For example, the following works without issues

```@repl rad; continued = true
s + s
```

```@repl rad; continued = true
v - 2*v
```

This next example, however, will fail since the shapes are incompatible:

```@repl rad; continued = true
show(size(s))
show(size(v))
s + v
```

When it so happens that two or more radiance objects have to be used in mathematical operations, a few simple checks can be done to make sure that compatible operations are performed.


```@repl rad; continued = true
new_rad = VectorRadiance(T, N); # T,N from above..
r1 = s; r2 = v; # or can use r1 = v; r2 = s;

@views new_rad[:,1] = r1[:,1] + r2[:,1]; # Add intensity component (both have them)

if r1 isa VectorRadiance
@views new_rad[:,2:end] += r1[:,2:end]
end

if r2 isa VectorRadiance
@views new_rad[:,2:end] += r2[:,2:end]
end
new_rad
```

In general, the expectation is that users tend to know which radiance objects are vector radiances, and which ones are scalar radiances, so it should generally be possible to write explicit and performant code when e.g. scalar solar (ir)radiances are multiplied with top-of-atmosphere vector radiances. For rare cases where the type is not known, type checks like above can be used to make sure the operations will succeed regardless of the type.