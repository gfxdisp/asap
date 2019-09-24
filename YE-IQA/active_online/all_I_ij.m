% This function evaluates I(\mathcal{E}_{ij},Pr(\theta)) for all i,j.
function IGs=all_I_ij(post_mean,post_cov,sigma)
    n=numel(post_mean);
    IGs=zeros(n);
    for i=1:n
        for j=1:(i-1)
            mu=[post_mean(i),post_mean(j)];
            omega=[post_cov(i,i),post_cov(i,j);post_cov(j,i),post_cov(j,j)];
            IGs(i,j)=I_ij(mu,omega,sigma);
        end
    end
end
    