# [Working with Dictionaries](@id working_with_dicts)

## Introduction

Dictionaries in Julia are very similar to dictionaries in other scripting languages, and most details can be found in the official documentation [here](https://docs.julialang.org/en/v1/base/collections/#Dictionaries). The basic idea is the same as in Python, for example: a dictionary is a collection of items that are accessed by *keys*. In general, almost anything that can be assigned as a variable in Julia is a valid item, and every hash-able object can act as a key.

One would write in Python
```python
# Python example
d = {1: "one", 2: "two"}
```

and similarly in Julia

```@example dict
# Julia example
d = Dict(1 => "one", 2 => "two")
```

In both Python and Julia, the values can be accessed in various ways, such as directly typing the key:

```python
# Python
d[1]
# should result in output of "one"
```

and here the Julia output:

```@example dict
# Julia
d[1]
```

There are differences between Python's and Julia's dictionaries however, and two of those are particularly important for RetrievalToolbox. First, many more of Julia's various objects are hash-able, and can therefore be used as keys in dictionaries. For example, we can create a dictionary whose keys are arrays:

```@example dict
d = Dict([1,2,3] => "first three", [4,5] => "four and five")
```

The example above is not possible in (pure) Python as lists or numpy arrays are not hash-able by default. The second important difference relates to the specific data type of dictionaries. In Julia, they are not simply "dictionaries", but can be *specific* dictionary types that lay out the specific types of both keys and items. In the example above, the dictionary type is actually printed along with the contents: this dictionary `d` maps a vector of 64-bit integers (`Vector{Int64}`) to strings (`String`).

Once created, it is possible to add new items to the collection **as long as they match the dictionary type**. Thus, we can add another item:

```@example dict
d[[100]] = "only one hundred"
d
```

However, the following fails, since the key is not a list of 64-bit integers:
```@repl dict
d[101] = "one hundred and one"
```

Julia attempts by default to automatically convert `101` into a `Vector{Int64}` object, but no such conversion exists, and therefore the call above fails.

In order to allow for key/item pairs of different types to be stored in one dictionary, one must initialize the dictionary accordingly. The most general form is a dictionary of type `Dict{Any,Any}` which can map any key to any item:

```@example dict
d = Dict{Any,Any}()
d["something"] = [1,2,3]
d['a'] = 1:3
d[123.456] = Float32
d
```

As can be seen in the code above, we can add key/item pairs of `String`/`Vector{Int64}`, `Char`/`UnitRange{Int64}` and `Float64`/`DataType` into the dictionary, and it will be accepted since `d` can use `Any` for both keys and items.

Usage of `Any` in dictionaries can be convenient, however is generally discouraged in cases where performance matters. In RetrievalToolbox, as we will see below, dictionaries are used for bookkeeping mostly, and there is no performance-critical component directly attached to them. There are a few cases in RetrievalToolbox, where `Any` is used as either the key or item type.

## Usage in RetrievalToolbox

Dictionaries in RetrievalToolbox are used to perform bookkeeping in order to connect different objects together that are related without having some explicit dependency. A very instructive example are objects of type `EarthAtmosphereOpticalProperties`. Inside those objects we find fields such as `aerosol_tau` or `gas_tau`. Those fields are intended to contain the spectrally dependent per-layer optical depths due to aerosols or gases present in the model atmosphere. The optical depth information is stored as array to allow for high-performance computations.

Bookkeeping is crucial; we want to minimize possible assignment errors, meaning the confusion of which gas and which array belong together. Traditionally, the bookkeeping would be done via some form of index which would be carried throughout the program to (for example) signify that index `1` represents the oxygen gas, and index `2` is water vapor.

Dictionaries allow us to make this bookkeeping more intuitive and more elegantly. the `gas_tau` field, for example, is a dictionary of type `Dict{GasAbsorber{T}, Array{T, 2}}`. In this dictionary we directly map a `GasObserver` object to the matching array. There is no need to carry an additional bookkeeping token which links the correct gas absorber to the correct array.

Another example is how radiative transfer (RT) objects are managed inside the `EarthAtmosphereBuffer`. Since it is possible to arrange for several spectral windows to be part of a buffer, the RT methods are indexed in a dictionary of type `Dict{<:AbstractSpectralWindow, <:AbstractRTMethod}`. For some `EarthAtmosphereBuffer` `buf`, and a spectral window `swin`, we can write

```julia
rt = buf.rt[swin]
```

which would return the `AbstractRTMethod` object that is logically linked to our spectral window. Note again here the distinction. The `AbstractRTMethod` object itself **does not have an explicit dependence** on the spectral window, even though its field `optical_properties` does. Hence, dictionaries like this are usually employed when the direct dependence is not established within object properties themselves, but needs to be tracked overall.

So internally, this type of referencing can reduce potential programming errors. There is also a benefit on the user side. If users want to inspect some underlying quantity, like the layer-resolved optical depths due to gas absorption, the also do not have to worry about some arbitrary index that links a gas to this quantity. One could just inquire via the gas objects themselves.

```julia
# let `swin` be the spectral window of interest, and `buf` be a buffer object
rt = buf.rt[swin] # Get the RT object from `swin`
gas_o2 = RE.get_gas_from_name(buf.scene.atmosphere, "O2") # Get the gas object
# We can now obtain the optical depth array for the gas `gas_o2`
gas_tau_o2 = rt.optical_properties.gas_tau[gas_o2]
```

In above code example, there was no need to use any arbitrary index to access the objects of interest, the objects themselves reference each other where appropriate: starting from some spectral window `swin`, we take the RT object `rt` via the dictionary `buf.rt`. Since we also have our gas observer `gas_o2`, we can directly obtain the gas optical depth array through the `rt.optical_properties` dictionary by use `gas_o2` as the *key*.

Note that adding allocated objects to an existing dictionary, or even creating a dictionary with already allocated objects, does not (in general) allocate new memory! Therefore, there is no additional memory footprint due to objects appearing in one or more dictionaries.

Look at the memory footprint below. The first call allocates a new array in memory and resevers several MiB of memory in the process. The second call only establishes the dictionary via referencing and does not need to allocate the existing array anew.

```@repl dict
_blank = Dict(1 => rand(2,2)) # hide
@time a = rand(500, 5_000);
@time d = Dict(1 => a);
```

!!! tip
    Creating a dictionary in Julia with existing objects generally does not allocate new memory!