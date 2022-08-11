"""
FreeWilsonSpectrum(N::Int64, Lam::Number, epsilon::Vector{<:Number}, t::Vector{<:Number})
single particle spectrum of the Wilson chain
without spins
# Output 
energy spectrum En in ascending order
"""
function FreeWilsonSpectrum(N::Int64, Lam::Number, epsilon::Vector{<:Number}, t::Vector{<:Number})
	if size(t,1) < N
		@error "Parameters smaller than iteration"
		return
	end
	H = zeros(N, N)
	H[1, 2] = t[1]
	H[1, 1] = epsilon[1]
	for k = 2:(N-1)
		H[k, k+1] = t[k]
		H[k, k-1] = t[k-1]
		H[k, k] = epsilon[k]
	end
	H[N, N-1] = t[N-1]
	H[N, N] = epsilon[N]
	
	Hnew = Lam^((N-1)/2)*H
	Es, Vs = eigen(Hnew)
	return sort(Es)
end

"""
zero field susceptibility from grand canonical ensemble
- 'N': number of Wilson chain site,
- 'Λ': logarithm discretization
- 'epsilon': onsite energy for Wilson chain electron site,
- 't': hopping strength for Wilson chain electron site,
- 'barbetaKWW': constant for temperature
# Output
Tχ temperature times magnetic susceptibility
T temperature points
"""
function ConductionChi(N::Int64, Λ::Number, epsilon::Vector{<:Number}, t::Vector{<:Number}; barbetaKWW::Number=0.6)
	Tχ = zeros(N-1)
	T = zeros(N-1) # temperature kB_T=1/beta
	barbeta = (1+1/Λ)/barbetaKWW/2
	for n = 2:N
		En = FreeWilsonSpectrum(n,Λ,epsilon,t)
		T[n-1] = Λ^(-(n-2)/2)*barbeta
		for k = 1:n
			Tχ[n-1] += 1/(8*cosh(En[k]*Λ^(-(n-1)/2)/(2*T[n-1])).^2)
		end
	end
	return Tχ, T
end

"""
Conduction band entropy
same as the conduction band chi
"""
function ConductionEnt(N::Int64, Λ::Number, epsilon::Vector{<:Number}, t::Vector{<:Number}; barbetaKWW::Number = 0.6)
        Sent = zeros(N-1)
        T = zeros(N-1) # temperature kB_T=1/beta
        barbeta = (1+1/Λ)/barbetaKWW/2
        for n = 2:N
                En = FreeWilsonSpectrum(n,Λ,epsilon,t)
                T[n-1] = Λ^(-(n-2)/2)*barbeta
                Fenergy, MeanE = 0, 0
		for k = 1:n
			betaE = En[k]*Λ^(-1/2)/barbeta
			if En[k] < -100
				Fenergy += 2*betaE
			else
				Fenergy -= 2*log(1+exp(-betaE))
			end
			MeanE += 2*betaE/(1+exp(betaE))
		end
		Sent[n-1] = MeanE - Fenergy
        end
        return Sent, T
end


"""
- 'Lambda': logarithm discretization
- 'GS':  energy eigenvalues from NRG,
- 'GSSz': Sz of MB eigenstate from NRG,
- 'barbetaKWW': constant for temperature
# Output
- 'Tchi' : magnetic susceptibility
- 'T' : temperature points
"""
function ChiImpFull(Λ::Number, GS::Vector{Vector{Float64}}, GSSz::Vector{Vector{Float64}};  barbetaKWW = 0.6)
	N = size(GS,1) + 1
	Sz2ave = zeros(N-1)
	Szave = zeros(N-1)
	ZN = zeros(N-1)
	T = zeros(N-1)
	barbeta = (1+1/Λ)/barbetaKWW/2
	for n = 2:N
		T[n-1] = Λ^(-(n-2)/2)*barbeta
		Ns = size(GS[n-1],1)
		for k = 1:Ns
			ZN[n-1] += exp(-GS[n-1][k]*Λ^(-(n-2)/2)/T[n-1])
			Sz2ave[n-1] += GSSz[n-1][k]^2*exp(-GS[n-1][k]*Λ^(-(n-2)/2)/T[n-1])
			Szave[n-1] += GSSz[n-1][k]*exp(-GS[n-1][k]*Λ^(-(n-2)/2)/T[n-1])
		end
		Sz2ave[n-1] /= ZN[n-1]
		Szave[n-1] /= ZN[n-1]
	end
	Chi = Sz2ave .- Szave.^2
	return Chi, T
end


"""
- 'Lambda': logarithm discretization
- 'GS':  energy eigenvalues from NRG,
- 'GSSz': Sz of MB eigenstate from NRG,
- 'barbetaKWW': constant for temperature
# Output
- 'Sent' : magnetic susceptibility
- 'T' : temperature points
"""
function SentropyFull(Λ::Number, GS::Vector{Vector{Float64}};  barbetaKWW = 0.6)
	N = size(GS,1)+1
	Eave = zeros(N-1)
	ZN = zeros(N-1)
	T = zeros(N-1)
	barbeta = (1+1/Λ)/barbetaKWW/2
	for n = 2:N
		T[n-1] = Λ^(-(n-2)/2)*barbeta
		Ns = size(GS[n-1],1)
		for k = 1:Ns
			Ek = GS[n-1][k]/barbeta
			ZN[n-1] += exp(-Ek)
			Eave[n-1] += Ek*exp(-Ek)
		end
		Eave[n-1] /= ZN[n-1]
	end
	Sent = Eave .+ log.(ZN)
	return Sent, T
end

	



