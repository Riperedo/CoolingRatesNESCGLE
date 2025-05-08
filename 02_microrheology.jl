using NESCGLE
using DelimitedFiles

# Include functions from the src directory
include("src/MicrorheologyFunctions.jl")

#############
#	main	#
#############

# cooling rates
CR = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001]

# directories
begin_folder = "Hysteresis/WCA/rate_"
end_folder = "/phi_0.700000/T_3.000000_to_0.000999"


# main loop
for α in CR
	# constructiong folders
	base_folder = begin_folder*complete_str(α)*end_folder
	ls_folder = filter(x -> startswith(x, "TP_"), readdir(base_folder))
	# reading waiting time-data
	time = readdlm(base_folder*"/waiting_times.dat", skipstart=1)[:,1]
	η = readdlm(base_folder*"/waiting_times.dat", skipstart=1)[:,5]
	idx = readdlm(base_folder*"/waiting_times.dat", skipstart=1)[:,end]
	temperature = infer_Temp(base_folder).(time)
	# reading each file into folders
	for filename in ls_folder
		#println(get_idx(filename), " ", get_idx_int(get_idx(filename)))
		tail = get_idx(filename)
		#idx = get_idx_int(tail)
		τ = readdlm(base_folder*"/"*filename, skipstart=1)[:,1]
		W = readdlm(base_folder*"/"*filename, skipstart=1)[:,6]
		z, G, G´, G´´ = Mason(τ, W)
		save_data(base_folder*"/micro_rheology_"*tail, [z G G´ G´´], header = "z G G´ G´´")
	end
	save_data(base_folder*"/hysteresis_process.dat", [time temperature η idx], header="time temperature η idx")
end
