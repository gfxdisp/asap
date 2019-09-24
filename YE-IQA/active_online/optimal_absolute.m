% Computes the optimal absolute judgement "suggestion" (to judge condition
% "suggestion") for a given multivariate normal posterior with mean
% post_mean and covariance post_cov. gamma_hat is the estimate for gamma.
% sigma denotes the standard deviation of performance. The function also
% optionally returns the expected information gain ig from the optimal
% judgement.
function [suggestion,ig]=optimal_absolute(post_mean,post_cov,gamma_hat,sigma)
    IGs=all_I_i(post_mean,post_cov,gamma_hat,sigma);
    [ig,suggestion]=max(IGs);
end