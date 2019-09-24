% This function evaluates I(\mathcal{E}_i,Pr(\theta)) for all i.
function IGs=all_I_i(post_mean,post_cov,gamma_hat,sigma)
    n=length(post_mean);
    IGs=zeros(1,n);
    for i=1:n
        mu0=post_mean(i);
        sigma0=post_cov(i,i);
        IGs(i)=I_i(mu0,sigma0,gamma_hat,sigma);
    end
end
    