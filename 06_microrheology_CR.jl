using NESCGLE
using DelimitedFiles

# Include functions from the src directory
include("src/MicrorheologyFunctions.jl")

#############
#	main	#
#############
const COOLING_RATES = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001]

# Base directory components for constructing paths
const FOLDER_PREFIX = "FiniteRate/WCA/rate_"
const FOLDER_SUFFIX = "/phi_0.700000/T_3.000000_to_0.000999" # This suffix forms part of the directory name

# main loop
for cooling_rate in COOLING_RATES
    # Construct the base directory path for the current cooling rate
    # complete_str(cooling_rate) is assumed to return the string representation of the rate for the path
    base_dir_path = FOLDER_PREFIX * complete_str(cooling_rate) * FOLDER_SUFFIX

    # List relevant files (e.g., "TP_...") in the base directory
    # Ensure base_dir_path actually exists before calling readdir if necessary
    tp_files = filter(x -> startswith(x, "TP_"), readdir(base_dir_path))

    # --- Process waiting_times.dat ---
    waiting_times_filepath = joinpath(base_dir_path, "waiting_times.dat")
    waiting_times_data = readdlm(waiting_times_filepath, skipstart=1)
    
    time_points = waiting_times_data[:,1]
    viscosity_η = waiting_times_data[:,5]
    indices_idx = waiting_times_data[:,end] # Assuming this is the correct 'idx'

    # Infer temperature based on time points and base directory
    temperatures = infer_Temp_CR(base_dir_path).(time_points)

    # --- Process each TP_ file ---
    for filename in tp_files
        file_path = joinpath(base_dir_path, filename)
        
        # Extract identifier from filename
        file_identifier_tail = get_idx(filename) # e.g., "001.dat" or "001"
        
        # Read data from the current TP_ file
        tp_data = readdlm(file_path, skipstart=1)
        τ_values = tp_data[:,1]
        W_values = tp_data[:,6]
        
        # Calculate Mason parameters
        z, G, G_prime, G_double_prime = Mason(τ_values, W_values)
        
        # Save microrheology results
        output_micro_rheology_path = joinpath(base_dir_path, "micro_rheology_" * file_identifier_tail)
        save_data(output_micro_rheology_path, [z G G_prime G_double_prime], header = "z G G_prime G_double_prime")
    end

    # Save aggregated cooling process data
    output_cooling_process_path = joinpath(base_dir_path, "cooling_rate_process.dat")
    save_data(output_cooling_process_path, [time_points temperatures viscosity_η indices_idx], header="time temperature η idx")
end
