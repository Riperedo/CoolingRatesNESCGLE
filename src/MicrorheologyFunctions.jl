using SpecialFunctions
using DelimitedFiles
"""
    complete_str(x)

Converts a number `x` to a string and pads it with trailing zeros
until it reaches a length of 8. If the string representation of `x`
is already 8 characters or longer, it truncates it to 8 characters.

This is used to create filenames consistent with a specific folder-format.

# Arguments
- `x`: The number to be converted and formatted.

# Returns
- `String`: A string of length 8.
"""
function complete_str(x)
    str = string(x)
    while length(str) < 8
        str *= "0"
    end
    return str[1:8]
end

"""
    der(τ, W, i)

Calculates the first and second derivatives of `log(W)` with respect to `log(τ)`
at point `i` using a three-point finite difference formula.

This is typically used for analyzing data where quantities are plotted on a log-log scale.

# Arguments
- `τ::AbstractVector`: Vector of time points (or a similar independent variable).
- `W::AbstractVector`: Vector of corresponding dependent variable values (e.g., Mean Squared Displacement).
- `i::Integer`: The index of the point in `τ` and `W` at which to calculate the derivatives.
                 `τ[i]` and `W[i]` are the central point, `τ[i-1], W[i-1]` and `τ[i+1], W[i+1]` are the neighbors.

# Returns
- `Tuple{Float64, Float64}`:
    - `m`: The first derivative `d(log(W))/d(log(τ))` at `τ[i]`.
    - `β`: The second derivative `d²(log(W))/d(log(τ))²` at `τ[i]`.
"""
function der(τ, W, i)
    f0 = log(W[i-1])
    f1 = log(W[i])
    f2 = log(W[i+1])
    x0 = log(τ[i-1])
    x1 = log(τ[i])
    x2 = log(τ[i+1])
    m = f0 * (x1 - x2)/((x0 - x1)*(x0 - x2)) + f1 * (2*x1 - x0 - x2)/((x1 - x0)*(x1 - x2)) + f2 * (x1 - x0)/((x2 - x0)*(x2 - x1)) 
    β = (2*f0)/((x1-x0)*(x2-x0))-(2*f1)/((x2-x1)*(x1-x0))+(2*f2)/((x2-x1)*(x2-x0))			
    return m, β
end

"""
    Mason(τ::AbstractVector, Δr²::AbstractVector)

Calculates viscoelastic moduli (G, G', G'') from Mean Squared Displacement (MSD) data
using the Mason-Weitz approximation.

The method relates the local power-law exponent of the MSD to the complex modulus.
The frequency `s` is taken as `1/τ`.

# Arguments
- `τ::AbstractVector{Float64}`: Time lags.
- `Δr²::AbstractVector{Float64}`: Mean Squared Displacement corresponding to `τ`.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}`:
    - `z`: Frequencies `s = 1/τ` (excluding points where derivative cannot be calculated).
    - `G`: Magnitude of the complex shear modulus `|G*(s)|` (approximate).
    - `G´`: Storage modulus `G'(s)`.
    - `G´´`: Loss modulus `G''(s)`.

The calculation uses:
- `α = d(log(Δr²))/d(log(τ))` (local logarithmic slope).
- `g = C / (Δr²(τ=1/s) * Γ(1+α(s)))`, where `C=6` is a constant prefactor.
- `G'(s) = g * cos(π*α(s)/2)`
- `G''(s) = g * sin(π*α(s)/2)`
"""
function Mason(τ, Δr²)
	N = length(τ)
	z = []
	G = []
	G´ = []
	G´´ = []
	for i in 3:N
		s = 1/τ[i]
		α = (log(Δr²[i-1])-log(Δr²[i]))/(log(τ[i-1])-log(τ[i]))
		g = 6/(Δr²[i]*gamma(1+α))
		g´ = g*cos(π*α/2)
		g´´ = g*sin(π*α/2)
		append!(z, s)
		append!(G, g)
		append!(G´, g´)
		append!(G´´, g´´)
	end
	return z, G, G´, G´´
end

"""
    get_idx(filename::String)

Extracts the last part of a filename, assuming it's an identifier separated by underscores.

# Arguments
- `filename::String`: The input filename.

# Returns
- `SubString{String}`: The substring after the last underscore.
                       For example, for "TP_001", it returns "001".
"""
function get_idx(filename::String)
	idx = split(filename, "_")
	return idx[end]
end

"""
    get_idx_int(filename::SubString)

Converts a filename (or a part of it, typically an index) to an integer.
It splits the input string by "." and parses the first part as an integer.

# Arguments
- `filename::SubString`: The input substring, usually an index extracted from a full filename
                         (e.g., "001.dat" or just "001").

# Returns
- `Int`: The parsed integer. For "001.dat", it returns 1. For "001", it returns 1.
"""
function get_idx_int(filename::SubString)
	idx = split(filename, ".")
	return parse(Int, idx[1])
end

"""
    infer_Temp(folder::String) -> Function

Infers a time-dependent temperature profile function `T(t)` based on a folder path string.
The folder path is expected to encode cooling rate (α), initial temperature (Ti),
and final temperature (Tf).

The temperature profile is assumed to be:
- `Ti - α*t` for `t <= tp` (cooling phase)
- `min(Tf + α*(t-tp), Ti)` for `t > tp` (potential reheating or stabilization phase, capped at Ti)
where `tp = (Ti-Tf)/α` is the time to reach Tf.

# Arguments
- `folder::String`: The folder path string. Expected format for relevant parts:
    - A subfolder like "rate_X" where X is the cooling rate `α`.
    - The last subfolder like "anything_Ti_to_Tf" where Ti is initial temperature
      and Tf is final temperature.
    Example: "Hysteresis/WCA/rate_0.1/phi_0.700000/T_3.0_to_0.001"

# Returns
- `Function`: A function `f(t::Float64)` that returns the temperature at a given time `t`.
"""
function infer_Temp(folder::String)
	subfolder = split(folder, "/")
	rate_folder = subfolder[3]
	Temp_folder = subfolder[end]
	α = parse(Float64, split(rate_folder, "_")[2])
	Ti = parse(Float64, split(Temp_folder, "_")[2])
	Tf = parse(Float64, split(Temp_folder, "_")[end])
	tp = (Ti-Tf)/α
	f(t::Float64) = t <= tp ? Ti - α*t : min(Tf + α*(t-tp), Ti)
	return f
end


function infer_Temp_CR(folder::String)
	subfolder = split(folder, "/")
	rate_folder = subfolder[3]
	Temp_folder = subfolder[end]
	α = parse(Float64, split(rate_folder, "_")[2])
	Ti = parse(Float64, split(Temp_folder, "_")[2])
	Tf = parse(Float64, split(Temp_folder, "_")[end])
	tp = (Ti-Tf)/α
	f(t::Float64) = t <= tp ? Ti - α*t : Tf
	return f
end

"""
    linear_interpolation(x_input::Vector{Float64}, y_input::Vector{Float64}, x::Float64)

Perform linear interpolation to estimate the value of `y` at a given point `x` based on input data points (`x_input`, `y_input`).

# Arguments
- `x_input`: A vector of x-coordinates of the input data points. Must be sorted in ascending order.
- `y_input`: A vector of y-coordinates corresponding to the input data points.
- `x`: The x-coordinate at which to estimate the value of `y`.

# Returns
The estimated value of `y` at the given `x` coordinate. If `x` is outside the range of `x_input`, the function returns the value of `y` at the nearest endpoint.

# Example
```julia
x_data = [1.0, 2.0, 3.0, 4.0]
y_data = [2.0, 4.0, 1.0, 3.0]
x_value = 2.5
y_estimated = linear_interpolation(x_data, y_data, x_value)  # Returns 2.5
```
"""
function linear_interpolation(x_input::Vector{Float64}, y_input::Vector{Float64}, x::Float64)
    if x <= x_input[1]
        printstyled("algandamanl\n", color=:red)
        return y_input[1]
    elseif x >= x_input[end]
        return y_input[end]
    else
        i = searchsortedfirst(x_input, x)
        x1, x2 = x_input[i-1], x_input[i]
        y1, y2 = y_input[i-1], y_input[i]
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    end
end

