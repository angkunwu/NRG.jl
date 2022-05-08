using NRG
using CSV
using DataFrames

# obtain biased mu through transformations
relmu = 0 #-0.7 # relative location for each transformation
relcenter = -0.16
mu0 = [relcenter]
trans = 10
for k=1:trans-1
	mu0[1] = mu0[1] +  relcenter/5^k	
end
μ = relmu/5.0^(trans)+mu0[1]
@show μ

Precision = 1000
setprecision(Precision)
Lam = BigFloat(5.0^(1/4))
Nmax = 4000
#μ = -0.80
Nwilson=40
maxlevel = 30
samescale = true
alphas, betas = NRG.ABwithmu(μ, Lam, Nmax; maxlevel,samescale)
@time t, epsilon = NRG.IntegralToWilsonParam(alphas,betas,Nwilson,Precision)


# export data using dataFrame to csv
df = DataFrame()
df.t = t
df.epsilon = epsilon
df.alphap = Float64.(alphas[1:Nwilson,1])
df.alphan = Float64.(alphas[1:Nwilson,2])
df.betap = Float64.(betas[1:Nwilson,1])
df.betan = Float64.(betas[1:Nwilson,2])

temp = "$μ" #[4:14]
temp = temp[4:end]
if ~samescale
	temp *="LRscales"
end
@show temp
filename = "KPMParamLam514mun0"*temp*"Nmax1k4Dig1kFractal.csv"
CSV.write(filename, df)

