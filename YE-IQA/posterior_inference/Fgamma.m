% Evaluates the function \mathcal{F}_\gamma(s) defined in the paper.
function F=Fgamma(s,mu,omega,gamma,P,M,sigma)
    F = -loglikeM(M, s, gamma, sigma)-loglikeP(P, s, sigma)-logprior(s,mu,omega);
end