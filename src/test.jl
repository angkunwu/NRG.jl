include("WilsonParam.jl")
include("FractalFillings.jl")


setprecision(300)
nmax = 100
N = 35
Lam = BigFloat(3.0)

alphas = zeros(BigFloat,nmax,2)
betas = zeros(BigFloat,nmax,2)

for k = 1:nmax
	alphas[k,1] = (Lam-1)/Lam^k/2
	betas[k,1] = alphas[k,1]*(Lam+1)/Lam^k/2
	alphas[k,2] = alphas[k,1]
	betas[k,2] = -betas[k,1]
end

t, epsilon = IntegralToWilsonParam(alphas, betas, N, 300);



Precision = 1000
setprecision(Precision)
Lam = BigFloat(5.0^(1/8))
Nmax = 4000
μ = 0.0
Nwilson=70
alphas, betas = ABwithmu(μ, Lam, Nmax)
@time t, epsilon = IntegralToWilsonParam(alphas,betas,Nwilson,Precision)



