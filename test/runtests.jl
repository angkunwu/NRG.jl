using Test, NRG
using DelimitedFiles

@show Threads.nthreads()

wilsons = readdlm("/Users/angkunwu/NRG/test/FlatWilsonParam.txt", ',', Float64)
ttest, epsilontest = wilsons[:,1], wilsons[:,2]

@testset "Wilson Params" begin
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
	Threads.@threads for k = 3:order
		mus[k] = (((-1)^(k-2)-1)/(k-2)-((-1)^k-1)/k)/4
	end
	@time alphas, betas = KPMmomentToIntegral(mus, Lam; Precision=Precision, nmax=nmax, z=z)

	alphas = zeros(BigFloat,nmax,2)
	betas = zeros(BigFloat,nmax,2)
	for k = 1:nmax
		alphas[k,1] = (Lam-1)/Lam^k/2
		betas[k,1] = alphas[k,1]*(Lam+1)/Lam^k/2
		alphas[k,2] = alphas[k,1]
		betas[k,2] = -betas[k,1]
	end
	@time t, epsilon = IntegralToWilsonParam(alphas, betas, N; Precision=Precision)

	@test t ≈ ttest
	@test epsilon ≈ epsilon
end
