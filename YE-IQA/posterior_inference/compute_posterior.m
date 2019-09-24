% Compute the mean and covariance matrix for the (multivariate normal)
% posterior Pr(s|P,M). 
% post_mean is \hat{s}_\hat{\gamma} in the paper. This is also the MAP
% estimate for s.
% post_cov is \hat{R}_\hat{\gamma}^{-1} in the paper.
function [post_mean,post_cov,gamma_hat] = compute_posterior(mu,omega,raw_M,sigma,raw_P,tau)
    % Add a small constant tau to the zero elements in M and P.
    P=aug_mat(raw_P,tau);
    M=aug_mat(raw_M,tau);
    
    % Compute the MLE for log(Pr(M,P|\gamma)).
    gamma_hat=compute_gamma_hat(P,M,sigma,mu,omega);
    
    % Compute the mean and covariance matrix for the (multivariate normal)
    % posterior.
    post_mean=compute_s_hat(mu,omega,gamma_hat,P,M,sigma);
    
    % This code can be used to normalize the posterior mean to a 0:(n-1)
    % scale.
    %post_mean=zscore(post_mean);
    %post_mean=post_mean-post_mean(1);
    %post_mean=post_mean/mean(post_mean(2:end)./(1:(length(post_mean)-1)));
    
    post_cov=inv(Rgamma(post_mean,mu,omega,gamma_hat,P,M,sigma));
end