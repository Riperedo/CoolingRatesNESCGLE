using DelimitedFiles
using Polynomials

#######################3

function save_data(nombre,formato)
	@assert typeof(nombre) == typeof("hola")
	open(nombre, "w") do io
		writedlm(io,formato)
	end
	println("Data saved as ", nombre)
end

function complete_str(x)
    str = string(x)
    while length(str) < 8
        str *= "0"
    end
    return str[1:8]
end

################################

##########################
##Modulos Viscoelasticos##
##########################
## G(ω) = -iωη(ω)

function Delhw2Gw(Omega, Delh_w)
	G1 = zeros(length(Omega))
	G2 = zeros(length(Omega))
	for i in 1:length(Omega)
		G1[i] = -Omega[i]*imag(Delh_w[i])
		G2[i] = Omega[i]*(real(Delh_w[i]) + 1)
	end
	return G1, G2
end

### TF de fourier para un polinomio de orden dos 
function FT_poly(A, B, C, a, b, k)
	return (exp(-1im*a*k)*(-1im*(k^2)*(a*(a*C + B) + A) - k*(2*a*C + B) + 2im*C) + exp(-1im*b*k)*(1im*(k^2)*(A + b*(b*C + B)) + k*(2*b*C + B) - 2im*C))/(k^3)
end

## Fit de orden dos (la emejor aproximacion es lineal por eso se utiliza como orden uno).
function polyfitL(x, y)
	P = fit(x, y, 1)
	A = P[0]
	B = P[1]
	C = P[2]
	return A, B, C
end

##TF fuerza bruta
function FTΔη(tau,δη,ω)
	dη = 0; N = length(tau)
	for i in 2:N
		dt = tau[i]-tau[i-1]
		dη += δη[i]*exp(-im*ω*tau[i])*dt 
	end
	return dη
end

##Fit de los Modulos Viscoelasticos Lineal
function SDFTL(Deta, tau, Omega)
	FT	= zeros(ComplexF64,length(Omega))	
	p = floor(Int64,length(Deta)/2)-1
	N = map(Int,length(Deta))
	for i in 0:p
		 if 2*i + 3 > N break end
		 	idx_a = 1 + 2*i
		 	idx_b = 3 + 2*i
			xs =	tau[idx_a:idx_b]
		 	ys = Deta[idx_a:idx_b]
			A, B, C= polyfitL(xs, ys)
			a = tau[idx_a]
			b = tau[idx_b]
			for i1 in 1:length(Omega)
					FT[i1] += FT_poly(A,B,C,a,b,Omega[i1])
			end
	end
	return FT
end


###########################333



function main(filename)
	global ηₜ = Float64[]
	global Deta = Float64[]
	global tol = 0; global steps = 1; global nsb = 3

	##Grid de frecuencias
	global dw = (10^-(1/5))
	global omega = [1e-12]
	for i in range(1,1,100)
	 	omega=append!(omega,omega[end]/dw)
	end 

	################################################
	##Viscosdiad como funcion del tiempo de espera##
	################################################
	data = readdlm(filename, skipstart = 1)
	τ = data[:,1]
	Δη = data[:,5]
	#visc = 1 .+ eta(τ,Δη) 
	#reduce(append!,(ηₜ, visc)) 
	#######################
	##Iniciando Reologia###
	#######################

	###########################
	##Modulos Visco-Elasticos##
	###########################
	BF = [FTΔη(τ[begin:steps:end], Δη[begin:steps:end], w) for w in omega]
	g_re, g_im = Delhw2Gw(omega,BF)
		
	i=1; flag=0
	while flag==0
		i=i+1
		dif1= g_re[i]-g_re[i-1]
		if	 dif1 < tol
		global ω_tol = omega[i-nsb]
		println(ω_tol)
		flag = 1
		end 
		end
		G_p_L = zeros(length(omega))
		G_b_L = zeros(length(omega))
		delh_w_L = SDFTL(Δη[begin:steps:end], τ[begin:steps:end], omega)
			G1, G2 = Delhw2Gw(omega, delh_w_L)
		for i2 in 1:length(omega)
		if omega[i2] < ω_tol
			G_p_L[i2] = g_re[i2]
			G_b_L[i2] = g_im[i2]
		 else
			G_p_L[i2] = G1[i2]
			G_b_L[i2] = G2[i2]
		end
	end
	return omega, G_p_L, G_b_L
end


#phi = LinRange(0.5, 0.6, 11)
phi = [0.5791, 0.5792, 0.5793, 0.5794, 0.5795, 0.5796, 0.5797, 0.5798, 0.5799]

for ϕ in phi
	filename = "dyn_eq_HS_phi_"*complete_str(ϕ)*".dat"

	omega, G_p_L, G_b_L = main(filename)

	save_data("rheology_true_eq_HS_phi_"*complete_str(ϕ)*".dat", [omega G_p_L G_b_L])
end