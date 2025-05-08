reset

set logscale y

set xlabel "T"
set ylabel "G''({/Symbol w})"

# Use the viridis color palette
set palette viridis

# Define constants for iteration
BASE_PATH_PREFIX = "FiniteRate/WCA/rate_"
BASE_PATH_SUFFIX = "/phi_0.700000/T_3.000000_to_0.000999/microrheology_cuts.dat"

# Cooling rates as they appear in directory names
RATE_STRINGS = "10.00000 1.000000 0.100000 0.010000 0.001000 0.000100"

# Corresponding exponents for the plot titles (10^exponent)
# 10.0    -> 10^1
# 1.0     -> 10^0
# 0.1     -> 10^-1
# 0.01    -> 10^-2
# 0.001   -> 10^-3
# 0.0001  -> 10^-4
EXPONENT_STRINGS = "1 0 -1 -2 -3 -4"

# Get the number of rates for calculating color fractions
N_RATES = words(RATE_STRINGS)

# Build the plot command dynamically
plot for [i=1:N_RATES] \
    sprintf("%s%s%s", BASE_PATH_PREFIX, word(RATE_STRINGS, i), BASE_PATH_SUFFIX) \
    using 2:4 with linespoints \
    pointtype 7 \
    linecolor palette frac (N_RATES > 1 ? (1.0 - (real(i)-1.0)/(real(N_RATES)-1.0)) : 0.5) \
    title sprintf("CR = 10^{%s}", word(EXPONENT_STRINGS, i))
