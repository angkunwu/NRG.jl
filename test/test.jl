#=
using DelimitedFiles
using NRG
using Plots
using LaTeXStrings

#include("NRGiterations.jl")

N = 30
Lam = 3.0
U = 0.1
ef = -U/2
V = 0.05
Ns = 4000
wilsons = readdlm("/Users/angkunwu/NRG/test/FlatWilsonParam.txt", ',', Float64)
t, epsilon = wilsons[:,1], wilsons[:,2]
h = 10^(-8)
GS, GSSz, Szimp = NRG.runNRG(N,Lam,epsilon,t,U,ef,V,Ns;h=h);
barbetaKWW = 0.6
ChiT, T = NRG.ChiImpFull(Lam,GS, GSSz; barbetaKWW)
ChiT0, T = NRG.ConductionChi(N, Lam, epsilon,t; barbetaKWW)
ChiImp = ChiT .- ChiT0

Sent, T = NRG.SentropyFull(Lam, GS)
S0, T = NRG.ConductionEnt(N, Lam, epsilon, t)
Simp = Sent .- S0

#plot(T, ChiImp, xaxis=:log, xlabel=L"T", ylabel=L"T\chi",legend=:bottomright)
#hline!([1/8 1/4], linestyle=:dash)

plot(T, Simp, xaxis=:log, xlabel=L"T", ylabel=L"S_{imp}",legend=:bottomright)
hline!([log(4) log(2)], linestyle=:dash)

ChiLoc, T = NRG.ChiLocal(Lam, GS, Szimp, h;Nrelax=100)
=#

using NRG

#include("WilsonParam.jl")
#include("FractalFillings.jl")

Precision = 300
setprecision(Precision)
nmax = 100
N = 35
Lam = BigFloat(3.0)
z = 0.0

# const Hybridization KPM moments
order = 1000
mus = zeros(order)
mus[1] = 1
for k = 3:order
	mus[k] = (((-1)^(k-2)-1)/(k-2)-((-1)^k-1)/k)/4
end
@time alphas, betas = NRG.KPMmomentToIntegral(mus, Lam, Precision, nmax, N, z)

alphas = zeros(BigFloat,nmax,2)
betas = zeros(BigFloat,nmax,2)
for k = 1:nmax
	alphas[k,1] = (Lam-1)/Lam^k/2
	betas[k,1] = alphas[k,1]*(Lam+1)/Lam^k/2
	alphas[k,2] = alphas[k,1]
	betas[k,2] = -betas[k,1]
end

@time t, epsilon = NRG.IntegralToWilsonParam(alphas, betas, N, 300)


#=
Precision = 1000
setprecision(Precision)
Lam = BigFloat(5.0^(1/8))
Nmax = 4000
μ = 0.0
Nwilson=70
alphas, betas = ABwithmu(μ, Lam, Nmax)
@time t, epsilon = IntegralToWilsonParam(alphas,betas,Nwilson,Precision)
=#


