% Evaluates log(Pr(M,P|\gamma)).
function loglike=loglikeMP(P,M,gamma,sigma,mu,omega)
    Pdims=size(P);
    n=Pdims(1);
    s_hat= compute_s_hat(mu,omega,gamma,P,M,sigma);
    R_hat=Rgamma(s_hat,mu,omega,gamma,P,M,sigma);
    % Note that for the next bit we must use logdet, not logabsdet.
    if n<100
        loglike=-Fgamma(s_hat,mu,omega,gamma,P,M,sigma)+(n/2)*log(2*pi)-0.5*log(det(R_hat));
    else
        loglike=-Fgamma(s_hat,mu,omega,gamma,P,M,sigma)+(n/2)*log(2*pi)-0.5*logdet(R_hat);
    end
end