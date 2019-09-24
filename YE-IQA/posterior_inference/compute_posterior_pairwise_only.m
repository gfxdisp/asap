% Compute the mean and covariance matrix for the (multivariate normal)
% posterior Pr(s|P) (when we have no absolute judgements). 
% post_mean is \hat{s}_\hat{\gamma} in the paper. This is also the MAP
% estimate for s.
% post_cov is \hat{R}_\hat{\gamma}^{-1} in the paper.
function [post_mean,post_cov] = compute_posterior_pairwise_only(mu,omega,sigma,raw_P,tau)
    n=size(raw_P,1);
    raw_M=zeros(4,n);
    [post_mean,post_cov] = compute_posterior(mu,omega,raw_M,sigma,raw_P,tau);
end
    