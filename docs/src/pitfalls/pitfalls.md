# Known pitfalls & issues

This section of the documentation deals with known pitfalls regarding the usage of the toolkit, as well as more structural issues. Due to the highly modular philosophy of this software, many critical portions of an algorithm implementation are up to the user and, by design, not part of the toolkit itself. Thus, it is easily possible to _break_ some of the intended program flow.

## Azimuthal angles and convention for use with XRTM



## Mutability of many RetrievalToolbox objects



## Angles and trigonometric functions with Unitful.jl

In Julia, calling trigonometric functions is straightforward, with angles being naturally considered to be in units of radiants:

```@repl
sin(pi/2)
```

There are also trigonometric functions that ingest angles in degrees, which are conveniently named the same but append a `d` at the end, e.g. `sin` becomes `sind`, or `cos` becomes `cosd`. The behavior is as expected

```@repl
sind(90.0)
```

When using the Unitful.jl package, one can attach appropriate units to angle values. The library also provides an implementation of the common trigonometric functions such that the appropriate unit conversions take place internally:
```@repl
using Unitful # hide

x = 90u"°"
sin(x)
cos(x)
```

Similarly, the `cosd` and `sind` functions as provided by Unitful.jl accept degree-valued quantities that provide the correct result. In the case of some degree-valued quantity, `x`, the result happens to the same:

```@repl
using Unitful # hide

x = 55u"°"
sin(x)
sind(x)

```

Looking closer into what functions are called, we observe that `sin` and `sind` provide the interface to the correct function when the argument is a degree-valued quantity. However if the argument is radians-valued, `sind` will not perform the appropriate conversion!

```@repl
using Unitful # hide

x = 45u"°"
y = (pi/4)u"rad"

code_lowered(sin, tuple(typeof(x)))
code_lowered(sind, tuple(typeof(x)))
code_lowered(sin, tuple(typeof(y)))
code_lowered(sind, tuple(typeof(y)))

```

The code block above shows that `sind`, when called with a radians-valued quantity will simply perform a `deg2rad` conversion and use the `sin` function on the argument. More importantly, no error will be thrown when such a computation is attempted!

```@repl
using Unitful # hide

y = (pi/4)u"rad"
z = y |> u"°"

# This is not what we want!
sind(y)
# This is!
sind(z)

```

!!! warning
    Trigonometric functions meant for degree-valued arguments, such as `sind` or `cosd` do not perform automatic unit conversions, even if the argument is a `Unitful` quantity! Calculating `cosd((pi/2)u"rad")` will yield an unexpected answer!


!!! tip
    It is best-practice to clearly document which units certain angle variables must have, and then use the appropriate trigonometric functions! Always use e.g. `cosd` when degree-valued angles are expected, but use e.g. `cos` if a `Unitful` angle quantity is expected.


Further, users should be aware of implicit conversions when assigning angle-valued quantities to object fields. This is a particular danger with buffers, as those are instantiated and then modified in-place at some later point. The example below illustrates the issue.

Assume some mutable, user-defined type `t` that accepts some float as its only field. We then instantiate a new variable `v` with some arbitrary value.

```@repl type; continued = true
using Unitful # hide
mutable struct t{T <: AbstractFloat}
    a::T
end

v = t(123.456)
```

At some later point, we would like to change `v.a`. For the sake of this example, assume we want to give it the value of 15 degrees. Following code executes without raising an error.

```@repl type

v.a = 15.0u"°"
```

Yet, when we inspect the value `v.a` itself, we see that an implicit conversion has taken place which turned the 15° into equivalent radians. The conversion happened since angular units are dimensionless, or have dimension [1], so can always be turned into a `Unitful.NoUnits` quantity which then can be cast into a "regular" Julia number. Some function which might perform a calculation on some object of this type, therefore should not use the degree-versions of trigonometric functions.

```@repl type
v.a
deg2rad(15.0)
```

!!! warning
    `Unitful` angle-valued quantities are prone to implicit conversions to radians!