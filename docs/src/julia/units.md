# Working with Unitful.jl units

## Introduction

The retrieval toolkit currently uses `Unitful.jl`[^Unitful] to attach physical units to various quantities that need them. This section here will present a short introduction to the package as well introduce helpful code snippets that should be used when working with RE. Some of these are currently used throughout the package to make sure that calculations are unit-aware.

`Unitful.jl` introduces new types which allows users to perform calculations that respect the units of quantities of those types. For example, adding various length-type quantities triggers automatic unit conversions (`using Unitful.jl` is assumed for all code examples).

```@repl
using Unitful # hide

5.2u"m" + 123.5u"inch" - 0.005u"mi"
```

Performing e.g. an addition with quantities of incompatible units will raise an error.

```@repl
using Unitful # hide
10.0u"m" + 5.0u"kg"
```

In general, most types of numbers and arrays can be outfitted with a physical unit, and most basic computations with them will involve the automatic checks for compatible dimensions and, if needed, unit conversions.

```@repl
using Unitful # hide
[1., 2., 3., 4.] * u"Pa" + [0.05, 0.06, 0.07, 0.08] * u"Torr"
```

Note that in the example above, the resulting quantity is not in units of Pascal, but in kg m⁻¹ s⁻², since `Unitful.jl` will usually default to SI units once a conversion has taken place. This behavior can be controlled, however, more details are found in the `Unitful.jl` documentation.

One can force unit conversion to any compatible unit via the `uconvert` function (note the broadcasting - `uconvert` does not dispatch on arrays or vectors).

```@repl
using Unitful # hide
x = [1., 2., 3., 4.] * u"Pa" + [0.05, 0.06, 0.07, 0.08] * u"Torr"
uconvert.(Ref(u"Pa"), x)
```

In the example above, it is important to understand that the object `x` is a vector of _pressure-type_ units. Julia allows for the mixing of types within a vector, so one can feasibly have a vector or array with mixed units. In this next example, we have three quantities of different units, although Julia will convert the quantities into `Float64`, even though the elements are initially typed as `Float64`, `Float32` and `Int64`.

```@repl
using Unitful # hide
[1.0u"m", 5.0f0u"kg", 20u"s"]
```

In general, computations with `Unitful`-type arrays tend to be almost as fast as those with "raw" arrays of type `Float32` or `Float64`. One of the exceptions is exponentiation with non-integer valued exponents, as shown below.

```@repl
using Unitful # hide
using BenchmarkTools

N = 16;
A = rand(Float32, N, N);
A_u = A * u"Pa";

@benchmark A.^3.5 seconds=1
@benchmark A_u.^3.5 seconds=1
```

Various linear algebra operations can be performed with `Unitful` quantities, such as the inverse of a square matrix, and the units of the returned object will be correct.

```@repl
using Unitful # hide

A = rand(5, 5) * u"kg";
inv(A)
```

On the other hand, computing the QR-decomposition on such a matrix will fail.

```@repl
using Unitful # hide
using LinearAlgebra
A = rand(5, 5) * u"kg";
qr(A)
```

!!! warning
    Many functions, especially those that call lower-level routines (such as BLAS), might not be compatible with `Unitful` quantities, or arrays composed of them!

While the result of certain linear algebra operations using matrices with physical units might make sense intuitively, others might not be well-defined. Further, it depends highly on the actual implementation and it is difficult to know whether a certain piece of Julia code will execute for a vector or matrix that has `Unitful` units attached to them.

## Usage within RetrievalToolbox

It is possible to use `Unitful` quantities directly in custom types and perform calculations using objects of those types. The convenience of this approach is clear: since every variable has units attached to them, an error will be raised immediately when an incompatible computation is performed. So if some part of a computation is accidentally mistyped and results in an operation that is not valid from a physical unit perspective, such as adding pressure and temperature, the code will not successfully execute. Further, if the computations themselves implicitly perform the required unit conversions, users do not need to write code that takes care of possible unit conversions. In addition to that, a flexible unit system allows users to choose their units of preference, e.g. state their pressures in "hPa" rather than "Pa".

There are currently two major drawbacks that were identified during the development of RetrievalToolbox when trying to use vectors and matrices of some `Unitful` type.

1. Some linear algebra operations were fundamentally incompatible with `Unitful` matrices
2. Some operations (e.g. exponentiation, taking the logarithm) suffered from performance drops

Regarding (1), most computations can be performed with `Unitful` vectors and matrices without having to adapt the code. Various linear algebra operations, however, as is necessary for the inversion, require matrices that by definition have mixed units. If `Unitful` is consistently used for all quantities, then the prior covariance matrix for an OE-type retrieval, for example, will have various units for the quantities that occur within the matrix. While the inversion equations obviously are consistent in terms of their physical units, some of the linear algebra operations (e.g. computing the inverse of a larger matrix might use decomposition techniques) will fail since those tend to be optimized for numerical efficiency and might contain incompatible sub-operations.
This issue can be circumvented by casting a `Unitful` matrix into a regular "numeric" matrix, performing the needed calculations, and then re-attaching the units at a later stage.

A similar solution can be used for (2) by identifying which operations suffer from performance issues and making sure that those will be executed for regular matrices and vectors, rather than involving the `Unitful` types.

The current implementation of physical units in RetrievalToolbox uses a hybrid and more explicit approach. Every custom type with a field that should have a physical unit, also has a field that contains the physical unit. As an example, a `GasAbsorber` contains a volume mixing ratio profile `vmr_levels`. We want to give users the option to define a VMR profile in various units, such as ppb when working with methane, or ppm when working with carbon dioxide. So for the `vmr_levels` field, there is an accompanying field `vmr_unit` which carries the `Unitful` type that describes the physical unit.

```julia
struct GasAbsorber{T} <: AbstractAtmosphereElement where T

    gas_name::String
    spectroscopy::AbstractSpectroscopy
    vmr_levels::Vector{T}
    vmr_unit::Unitful.DimensionlessUnits

end
```

Note that above, we can not just demand that the user instantiate the object with a `Unitful` type, but we can also make sure that the assigned unit makes sense for this particular object. Volume mixing ratios should be dimensionless (or dimension 1), and `Unitful` provides a type that represents all the explicit types. There are many other types that fall into the same category, we can use e.g. `Unitful.LengthUnits` if we want to force any type of length unit.

Of course now, users and code maintainers must use explicitly use the unit information and make sure that possible unit conversions are taken into account when they need to be. For dimensionless units, the situation is straightforward, since applying the `NoUnits` function will convert the dimensionless quantity into parts.

Users and code maintainers must ensure the appropriate use, and below we demonstrate the general idea with examples.

For example one, we have two dimensionless quantities `x` and `y`. Both `x` and `y` are `Float64` in this case, and their associated units are stored in separate variables `x_u` and `y_u`. If we have to compute a new quantity `x + y`, we must first make sure the units of each are considered.

```@repl
using Unitful # hide

x = 1.234;
x_u = u"ppm";

y = 1500.0;
y_u = u"ppb";

y * y_u + x * x_u |> NoUnits
```

In the next example, we consider two pressure values `p1` and `p2` whose units `p1_u` and `p2_u` are not necessarily the same. Let us assume that we must carry on the computation in units of `p2_u`. The function `ustrip` provided by `Unitful` is very helpful in this way. When called with a `Unitful` quantity alone, the unit will be simply stripped, and the numeric value is retained, e.g. `ustrip(123.0u"hPa")` returns `123.0`. However, when we call `ustrip` with a leading argument, being a compatible unit, the quantity will first be converted to that unit. `ustrip(u"Pa", 123.0u"hPa")` will return `12300.0`.

In the first variant of performing `p1 + p2`, we first construct `Unitful` quantities by multiplying each pressure value by its respective unit, `p1 * p1_u` and `p2 * p2_u`, and then perform the sum. The intermediate sum will be in SI pressure units, but then we apply the `ustrip` function to convert the sum into units of `p2_u` and strip the unit to obtain the numerical value only.

```@repl
using Unitful # hide
#=
    Version 1
=#
p1 = 1100.0; p1_u = u"hPa";
p2 = 55.5; p2_u = u"Pa";
result = ustrip(p2_u, p1 * p1_u + p2 * p2_u)
```

In the second variant demonstrated here we first construct `p1 * p1_u` to obtain a `Unitful` object, which yields 1100.0 hPa, and immediately convert it into units of `p2_u` (Pa) and stripping the unit, such that this new variable `tmp` becomes `110000.0`. Now we have pressure `p1` as a numerical value in units of Pa (implicit), and can simply add it to `p2` to obtain the same result as above.

```@repl
using Unitful # hide
#=
    Version 2
=#
p1 = 1100.0; p1_u = u"hPa";
p2 = 55.5; p2_u = u"Pa";
tmp = ustrip(p2_u, p1 * p1_u);
result = p2 + tmp
```

Throughout RetrievalToolbox both variants are applied.

[^Unitful]: [https://github.com/PainterQubits/Unitful.jl](https://github.com/PainterQubits/Unitful.jl)