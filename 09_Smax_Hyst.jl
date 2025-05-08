using NESCGLE
using DelimitedFiles # Es buena práctica incluirlo si usas readdlm explícitamente


# Include functions from the src directory
include("src/MicrorheologyFunctions.jl")

# Cooling rates (usar const para variables globales que no cambian)
const CR = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001]

# Directory parts (usar const)
const BEGIN_FOLDER = "Hysteresis/WCA/rate_"
const END_FOLDER = "/phi_0.700000/T_3.000000_to_0.000999"

# Main loop
for α in CR
	# Constructing folders
	base_folder = BEGIN_FOLDER * complete_str(α) * END_FOLDER
	ls_folder = filter(x -> startswith(x, "SF_"), readdir(base_folder))
	
	# Reading waiting time-data once
    waiting_data_path = joinpath(base_folder, "waiting_times.dat")
    waiting_data = readdlm(waiting_data_path, skipstart=1)
	time_points = waiting_data[:,1]
	original_indices = waiting_data[:,end] # Estos son los índices de archivo originales (probablemente Float64)
	
	temperature = infer_Temp(base_folder).(time_points)

    S_max = similar(time_points, Float64)

	# Reading each microrheology data file
	for filename in ls_folder
		tail = get_idx(filename)
		current_file_idx_int = get_idx_int(tail) # Índice entero del archivo actual

        # Encontrar la fila correspondiente en time_points, temperature, etc.
        # Es importante que la comparación sea robusta (ej. Int vs Float64)
        target_row_idx = findfirst(x -> isequal(x, Float64(current_file_idx_int)), original_indices)

        if target_row_idx === nothing
            @warn "Índice $current_file_idx_int del archivo $filename no encontrado en $waiting_data_path para la carpeta $base_folder. Omitiendo."
            continue
        end

        # Leer datos del archivo actual
        filepath = joinpath(base_folder, filename)
		file_data = readdlm(filepath, skipstart=1)
		S_data = file_data[:,2]

        S_max[target_row_idx] = maximum(S_data)

    end

    output_filename = joinpath(base_folder, "S_max.dat")

	save_data(output_filename, [time_points temperature S_max], header="t temperature Smax")
    println("Procesado y guardado datos para α = $α en $output_filename")
end

println("Todo el procesamiento completado.")
