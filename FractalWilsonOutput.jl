using NRG
using CSV
using DataFrames

# get μ through center indexes
inds = [1, 3, 2] # 1, left band; 2, middle band; 3 right band
length = size(inds,1)
bandleft, bandright = [-1.0], [1.0]
for k in 1:length
	bandwidth = bandright[1]-bandleft[1]
	center = (bandleft[1]+bandright[1])/2.0
	if inds[k] == 1
		temp = [center-0.5*bandwidth, center-0.3*bandwidth]
	elseif inds[k] == 2
		temp = [center-0.1*bandwidth, center+0.1*bandwidth]
	else
		temp = [center+0.3*bandwidth, center+0.5*bandwidth]
	end
	bandleft[1] = temp[1]
	bandright[1] = temp[2]
end
bandwidth = bandright[1]-bandleft[1]
center = (bandleft[1]+bandright[1])/2.0
μ = center + 0.11 * bandwidth 
@show μ


Precision = 1000
setprecision(Precision)
Lam = BigFloat(3.0) #BigFloat(5.0^(1/4))
Nmax = 4000
#μ = -0.16
Nwilson=40
maxlevel = 3
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

if μ ≠ 0
	temp = "$μ" #[4:14]
	temp = temp[4:end]
else
	temp = "0"
end
if ~samescale
	temp *="LRscales"
end
temp = temp*"mlev"*"$maxlevel"
@show temp
#filename = "KPMParamLam514mun0"*temp*"Nmax1k4Dig1kFractal.csv"
filename = "testLam3mun"*temp*".csv"
CSV.write(filename, df)

