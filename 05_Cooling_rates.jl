using NESCGLE

# --- System Parameters ---
const Nk = 200
const kmax = 15*π
const dk = kmax/Nk # Derived constant

# Wavevector k
const k = dk .* collect((1:Nk) .- 0.5)

# --- Initial and Final State Parameters ---
const ϕi = 0.7 # Initial packing fraction
const ϕf = 0.7 # Final packing fraction (same as initial in this case)
const Ti = 3.0 # Initial temperature
const Tf = 0.001 # Final temperature

const initial_SM = SM_WCA(ϕi, Ti, k)
const final_SM = SM_WCA(ϕf, Tf, k)

# --- Cooling Protocol Parameters ---
const COOLING_RATES = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001] # Renamed from CR
const NUM_PROTOCOL_STEPS = 40 # Number of steps for the FiniteRate protocol

for cooling_rate in COOLING_RATES # Renamed from α
    preparation_protocol = FiniteRate(initial_SM.params, final_SM.params, cooling_rate, NUM_PROTOCOL_STEPS)
    solution = NESCGLEsolver(initial_SM, preparation_protocol)
    save_files(solution, initial_SM, preparation_protocol)
end
