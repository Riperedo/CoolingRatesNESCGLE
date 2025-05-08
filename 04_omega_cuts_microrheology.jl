using NESCGLE
using DelimitedFiles # Es buena práctica incluirlo si usas readdlm explícitamente

# Include functions from the src directory
include("src/MicrorheologyFunctions.jl")

# Cooling rates (usar const para variables globales que no cambian)
const CR = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001]

# Omega cut values for interpolation
const OMEGA_CUTS = [1.0, 0.1, 0.01, 0.001, 0.0001]
const NUM_OMEGA_CUTS = length(OMEGA_CUTS)

# Directory parts (usar const)
const BEGIN_FOLDER = "Hysteresis/WCA/rate_"
const END_FOLDER = "/phi_0.700000/T_3.000000_to_0.000999"

# Main loop
for α in CR
	# Constructing folders
	base_folder = BEGIN_FOLDER * complete_str(α) * END_FOLDER
	ls_folder = filter(x -> startswith(x, "micro_rheology_"), readdir(base_folder))
	
	# Reading waiting time-data once
    waiting_data_path = joinpath(base_folder, "waiting_times.dat")
    waiting_data = readdlm(waiting_data_path, skipstart=1)
	time_points = waiting_data[:,1]
	original_indices = waiting_data[:,end] # Estos son los índices de archivo originales (probablemente Float64)
	
	temperature = infer_Temp(base_folder).(time_points)

    # Inicializar arrays para G' y G'' de forma programática
    # G_prime_at_cuts[i] almacenará los valores de G' para OMEGA_CUTS[i]
    G_prime_at_cuts = [similar(time_points, Float64) for _ in 1:NUM_OMEGA_CUTS]
    G_double_prime_at_cuts = [similar(time_points, Float64) for _ in 1:NUM_OMEGA_CUTS]

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
		ω_data = reverse(file_data[:,1])
		G_prime_data = reverse(file_data[:,3])
        G_double_prime_data = reverse(file_data[:,4])

        # Interpolar G' y G'' para cada valor de omega_cut y almacenar
        for k in 1:NUM_OMEGA_CUTS
            omega_val = OMEGA_CUTS[k]
            G_prime_at_cuts[k][target_row_idx] = linear_interpolation(ω_data, G_prime_data, omega_val)
            G_double_prime_at_cuts[k][target_row_idx] = linear_interpolation(ω_data, G_double_prime_data, omega_val)
        end
    end

    # Preparar datos y encabezado para guardar
    output_data_matrix = [time_points temperature] # Empezar con tiempo y temperatura
    header_parts = ["time", "temperature"]

    for k in 1:NUM_OMEGA_CUTS
        output_data_matrix = hcat(output_data_matrix, G_prime_at_cuts[k], G_double_prime_at_cuts[k])
        push!(header_parts, "G'($(OMEGA_CUTS[k]))", "G''($(OMEGA_CUTS[k]))")
    end
    
    full_header = join(header_parts, " ")
    output_filename = joinpath(base_folder, "microrheology_cuts.dat")

	save_data(output_filename, output_data_matrix, header=full_header)
    println("Procesado y guardado datos para α = $α en $output_filename")
end

println("Todo el procesamiento completado.")
