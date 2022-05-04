using NRG
using CSV
using DataFrames

# obtain biased mu through transformations
relmu = -0.7 # relative location for each transformation
relcenter = -0.8
mu0 = [-0.8]
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
#μ = -0.8
Nwilson=70
alphas, betas = NRG.ABwithmu(μ, Lam, Nmax)
@time t, epsilon = NRG.IntegralToWilsonParam(alphas,betas,Nwilson,Precision)


# export data using dataFrame to csv
df = DataFrame()
df.t = t
df.epsilon = epsilon
temp = "$μ"
filename = "KPMParamLam514mun"*temp[2:14]*"Nmax1k4Dig1kFractal.csv"
CSV.write(filename, df)

