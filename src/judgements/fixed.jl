"""
Fixed{T <: Number}

A wrapper type for fixed parameters in a measurement model. The `Fixed` type is used to indicate that a certain numeric parameter, such as a factor variance or loading in a `Measurement`, is fixed. 

# Arguments
- `x::T`: The numeric value being wrapped as fixed.

# Usage
Use to indicate fixed parameters in a measurement model.
To create a fixed parameter, simply wrap the numeric value with `Fixed`. For example, `Fixed(1.0)` creates a fixed value of `1.0`.
The strip_fixed function is provided to extract the numeric value from a Fixed object. It can also handle regular numeric types, returning the value as is.

```jldoctest
fixed_loading = Fixed(0.7)

# output

Fixed{Float64}(0.7)
```
"""
struct Fixed{T <: Number} <: Number
    value::T
end
strip_fixed(x::Fixed) = x.value
strip_fixed(x::Number) = x