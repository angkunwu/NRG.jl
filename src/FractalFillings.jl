function findinteraction(vec1, vec2)
    """
    given two ranges, return the interaction of the two range 
    """
    if vec1[2] < vec2[1] || vec2[2] < vec1[1]
        return []
    else
        return [max(vec1[1],vec2[1]), min(vec1[2],vec2[2])]
    end
end

function obtainAlphaBeta(rangevec, fractalband, level; maxlevel = 30)
    """
    given range [left, right] and fractalband for specific level
    return the estimated alpha, beta integrals falling in the fractalband
    up to l=30 levels
    """
    height = 0.5*(5.0/3)^level
    temp = findinteraction(rangevec, fractalband)
    if temp == fractalband
        alpha = 1/3^(level)
        beta = 0.5 * 2/5^level * sum(fractalband) * height
        return alpha, beta
    end

    if level == maxlevel
        alpha = 1/3^(level)
        beta = 0.5 * 2/5^level * sum(fractalband) * height
        return alpha, beta
    end

    level += 1
    width = 2/5^level
    rangel = [fractalband[1], fractalband[1]+width]
    rangec = [fractalband[1]+2*width, fractalband[2]-2*width]
    ranger = [fractalband[2]-width, fractalband[2]]

    alpha = 0
    beta = 0
    temp = findinteraction(rangevec, rangel)
    if ~isempty(temp)
        nextalpha, nextbeta = obtainAlphaBeta(temp, rangel, level)
        alpha += nextalpha
        beta += nextbeta
    end
    temp = findinteraction(rangevec, rangec)
    if ~isempty(temp)
        nextalpha, nextbeta = obtainAlphaBeta(temp, rangec, level)
        alpha += nextalpha
        beta += nextbeta
    end
    temp = findinteraction(rangevec, ranger)
    if ~isempty(temp)
        nextalpha, nextbeta = obtainAlphaBeta(temp, ranger, level)
        alpha += nextalpha
        beta += nextbeta
    end

    return alpha, beta
end


function ABwithmu(μ, Lam, N)
    """
    compute alphas betas with some chemical potential
    return alphas betas as big float
    """
    alphas = zeros(BigFloat, N, 2)
    betas = zeros(BigFloat, N, 2)
    fractalband=[-1, 1]
    level = 0
    scale = (1+abs(μ))
    for k in 1:N
        rangevec = [scale*Lam^(-k)+μ, scale*Lam^(1-k)+μ]
        rangevecn = [μ-scale*Lam^(1-k), μ-scale*Lam^(-k)]
        alphas[k, 1], betas[k, 1] = obtainAlphaBeta(rangevec, fractalband, level)
        alphas[k, 2], betas[k, 2] = obtainAlphaBeta(rangevecn, fractalband, level)
    end
    return alphas, betas
end 


