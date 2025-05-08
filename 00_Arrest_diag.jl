using NESCGLE
#using JSON

function main(args...)
	# user inputs
	ν_string = args[1]

	ν_string = split(ν_string, ".")[1]

	# variable definition
	ν = parse(Int64, ν_string)
            
    # preparing ϕ-T space grid
	phi = collect(0.581:0.001:0.8)
	Temperature = zeros(length(phi))
	T_min = 1e-10
	T_max = 1e1

	# preparing Input objects
	Nk = 200; kmax = 15*π; dk = kmax/Nk
	k = dk*(collect(1:Nk).- 0.5)

	# main loop
	for (i, ϕ) in enumerate(phi)

		function condition(T)
			sm = SM_WCA(ϕ, T, k, ν = ν)
			S = structure_factor(sm)
			iterations, gammas, state = Asymptotic(ϕ, k, S)
			return state == "Glass"
		end

		Temperature[i] = NESCGLE.bisection(condition, T_min, T_max, 1e-8)
	end
    save_data("Arrest_diag_wca.dat", [phi Temperature], header="phi temperature")
end

main(ARGS...)
# julia Arrest_diag_wca.jl 6