% Evaluates log(Pr(M|s,\gamma)).
% Note that the paper uses the notation M_{i,j} to denote the j^{th}
% row, i^{th} column.
function loglike = loglikeM(M, s, gamma, sigma)
    invsqrtsigma=sigma^(-0.5);
    gamma=gamma*invsqrtsigma;
    s=s*invsqrtsigma;
    
    diffmat=gamma.'-s;
    diffmat2=normcdf(diffmat(2:end,:))-normcdf(diffmat(1:(end-1),:));
    diffmat2(diffmat2==0)=1e-300;
    loglike=sum(sum(M.*log(diffmat2)));
end
