using NESCGLE

Nk = 200; kmax = 15*π; dk = kmax/Nk
k = dk*(collect(1:Nk) .- 0.5)
ϕi = 0.7
ϕf = 0.7
Ti = 3.0
Tf = 0.001
sm = SM_WCA(ϕi, Ti, k)
ℇ = SM_WCA(ϕf, Tf, k)
CR = [10.0, 1.0, 0.1, 0.01, 0.001, 0.0001]
for α in CR
    pp = Hysteresis(sm.params, ℇ.params, α, 30)
    sol = NESCGLEsolver(sm, pp)
    save_files(sol, sm, pp)
end
