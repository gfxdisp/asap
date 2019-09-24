% Compute the mean and covariance matrix for the (multivariate normal)
% posterior Pr(s|M) (when we have no pairwise judgements). 
% post_mean is \hat{s}_\hat{\gamma} in the paper. This is also the MAP
% estimate for s.
% post_cov is \hat{R}_\hat{\gamma}^{-1} in the paper.
function [post_mean,post_cov,gamma_hat] = compute_posterior_absolute_only(mu,omega,raw_M,sigma,tau)
    n=size(raw_M,2);
    raw_P=zeros(n);
    [post_mean,post_cov,gamma_hat] = compute_posterior(mu,omega,raw_M,sigma,raw_P,tau);
end
    