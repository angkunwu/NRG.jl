"""
        Use KPM moment to compute logarithmic integrals with Jackson Kernel
        With build-in paralization
        Precision: precision of numerics
        Lam: Logarithmic discretization parameter
        nmax: Set max # of logarithmic intervals or min range
        N: maximum Wilson chain sites
        z: z average z in [0,1)
"""
function KPMmomentToIntegral(mu::Vector{<:Number}, Lam::Number, Precision::Int64, nmax::Int64, N::Int64, z::Number)
    	order = size(mu, 1) 
    	BinaryPrec = Int(round(log(10)*Precision/log(2))) # convert precision to binary precision
    	setprecision(BinaryPrec)

    	g = zeros(BigFloat, order);	# Jackon kernel
    	for k = 1:order
		g[k] = ((order-k+2)*cos(pi*(k-1)/(order+1))+sin(pi*(k-1)/(order+1))*cot(pi/(order+1)))/(order+1)
    	end
	
    	epsilon = zeros(BigFloat,N); t = zeros(BigFloat,N);
    	#xi0=vpa(g(1)*mu(1));  #xi0=xi0.*(xi0>Threshold);
    	#epsilon(1)=vpa(g(2)*mu(2)); # Ekin at first site

    	theta = zeros(BigFloat,nmax+1);
    	theta[1]=acos(1);
    	for k=2:nmax+1
        	theta[k]=acos(Lam^(z-k+1));
    	end

	alphas = zeros(BigFloat, nmax, 2)
	betas = zeros(BigFloat, nmax, 2)
	# need paralization
    	#parfor l=1:nmax
    	for l=1:nmax
        	g21 = g[1]*mu[1]*(theta[l+1]-theta[l])
        	Threads.@threads for k=2:order
            		g21 += 2*g[k]*mu[k]*(sin((k-1)*theta[l+1])-sin((k-1)*theta[l]))/(k-1)
        	end
        	g21 /= pi
		alphas[l, 1] = g21
    
        	g22 = g[1]*mu[1]*(theta[l+1]-theta[l])
        	Threads.@threads for k = 2:order
            		g22 += 2*(-1)^(k-1)*g[k]*mu[k]*(sin((k-1)*theta[l+1])-sin((k-1)*theta[l]))/(k-1)
        	end
        	g22 /= pi
		alphas[l, 2] = g22

        	beta1 = g[2]*mu[2]*(theta[l+1]-theta[l])
        	Threads.@threads for k = 1:order-2
            		beta1 += (g[k]*mu[k]+g[k+2]*mu[k+2])*(sin(k*theta[l+1])-sin(k*theta[l]))/k
        	end
        	beta1 += g[order-1]*mu[order-1]*(sin((order-1)*theta[l+1])-sin((order-1)*theta[l]))/(order-1)
        	beta1 += g[order]*mu[order]*(sin(order*theta[l+1])-sin(order*theta[l]))/order
		betas[l, 1] = beta1./pi
    
        	beta2=g[2]*mu[2]*(theta[l+1]-theta[l]);
        	Threads.@threads for k=1:order-2
            	beta2 += (-1)^(k)*(g[k]*mu[k]+g[k+2]*mu[k+2])*(sin(k*theta[l+1])-sin(k*theta[l]))/k;
        	end
        	beta2 += (-1)^(order-1)*g[order-1]*mu[order-1]*(sin((order-1)*theta[l+1])-sin((order-1)*theta[l]))/(order-1);
        	beta2 += (-1)^(order)*g[order]*mu[order]*(sin(order*theta[l+1])-sin(order*theta[l]))/order;
		betas[l, 2] = beta2./pi
    	end

	return alphas, betas
end

"""
        From logarithmic grid integrals to the Wilson chian 
        
        alphas are the α integrals (2 vectors) with nmax rows
        betas are the β integrals (2 vectors) with nmax rows
        N is the number of the Wilson chain sites
        Precision is the high precision arithmetics
"""
function IntegralToWilsonParam(alphas::Matrix{<:BigFloat},betas::Matrix{<:BigFloat},N::Int,Precision::Int)
	if size(alphas,2) > 2 || size(betas,2) > 2 || size(alphas) != size(betas)
		error("αs or βs dimensions error!")
	end
	nmax = size(alphas,1) # number of logarithmic grid
	BinaryPrec = Int(round(log(10)*Precision/log(2))) # convert precision to binary precision
	setprecision(BinaryPrec)

	gamma2 = alphas
	xi = zeros(BigFloat, nmax, 2)
	for k = 1:nmax
		if alphas[k, 1] == 0
			xi[k, 1] = 0
		else
			xi[k, 1] = betas[k, 1]./ alphas[k, 1]
		end
		if alphas[k, 2] == 0
			xi[k, 2] = 0
		else
			xi[k, 2] = betas[k, 2]./ alphas[k, 2]
		end 
	end
	xi0 = sum(alphas)

	t = zeros(BigFloat, N)
        epsilon = zeros(BigFloat, N)
	
	u = zeros(BigFloat,N+1,nmax)
	v = zeros(BigFloat,N+1,nmax)
	
	epsilon[1] = sum(betas)/xi0 # xi0 == sum(alphas) == 1
	t[1]=sqrt((((xi[:,1].-epsilon[1]).^2)'*gamma2[:,1]+((xi[:,2].-epsilon[1]).^2)'*gamma2[:,2])/xi0);
	
	u[1,:] = sqrt.(gamma2[:,1]./xi0)';
	v[1,:] = sqrt.(gamma2[:,2]./xi0)';
	u[2,:] = (xi[:,1].-epsilon[1]).*u[1,:]./t[1];
	v[2,:] = (xi[:,2].-epsilon[1]).*v[1,:]./t[1];
	
	# Iterative solve for more parameters
	for k=2:N
		epsilon[k] = xi[:,1]'*(u[k,:]).^2+xi[:,2]'*(v[k,:]).^2;
		temp = BigFloat(0);
		for m=1:nmax
			temp += ((xi[m,1]-epsilon[k])*u[k,m]-t[k-1]*u[k-1,m])^2+((xi[m,2]-epsilon[k])*v[k,m]-t[k-1]*v[k-1,m])^2;
		end
		t[k] = sqrt(temp);
		u[k+1,:] = ((xi[:,1].-epsilon[k]).*u[k,:]-t[k-1]*u[k-1,:])./t[k];
		v[k+1,:] = ((xi[:,2].-epsilon[k]).*v[k,:]-t[k-1]*v[k-1,:])./t[k];
	end
	t = Float64.(t);
	epsilon = Float64.(epsilon);
	for k=1:N
		if abs(epsilon[k]) < 1E-100
			epsilon[k] = 0.0
		end
	end
	return (t,epsilon)
end










