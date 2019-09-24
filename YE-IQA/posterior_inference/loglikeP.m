% Evaluates log(Pr(P|s,\gamma)).
function loglike = loglikeP(P, s, sigma)
    s=s/sqrt(sigma*sqrt(2));
    cdfmat=normcdf(s.'-s);
    cdfmat(cdfmat==0)=1e-300;
    loglike=sum(sum(P.*log(cdfmat)));
end
