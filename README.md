# NRG.jl

|**Citation**                                                                     |**Open-access preprint**                               |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------:|
| | [![arXiv](https://img.shields.io/badge/arXiv-2205.06264-b31b1b.svg)](https://arxiv.org/abs/2205.06264) |


A julia package in numerical renormalization group with KPM integration

## Installation

The package is currently unregistered. To install it, first install and start Julia, then run the following command:
```julia
julia>]

pkg> add https://github.com/angkunwu/NRG.jl
```

Currently, you need to have permission of this package and enter Github username and PAT.

## Examples

Your can quickly test the KPM moment to the Wilson chain parameters for the matelic case
```julia
using NRG

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
@time alphas, betas = KPMmomentToIntegral(mus, Lam; Precision=Precision, nmax=nmax, z=z)

alphas = zeros(BigFloat,nmax,2)
betas = zeros(BigFloat,nmax,2)
for k = 1:nmax
	alphas[k,1] = (Lam-1)/Lam^k/2
	betas[k,1] = alphas[k,1]*(Lam+1)/Lam^k/2
	alphas[k,2] = alphas[k,1]
	betas[k,2] = -betas[k,1]
end

@time t, epsilon =  IntegralToWilsonParam(alphas, betas, N; Precision=Precision)
```


