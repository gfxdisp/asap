% Computes the optimal pairwise comparison suggestion=[a,b] (to compare a
% with b) for a given multivariate normal posterior with mean post_mean
% and covariance post_cov. sigma denotes the standard deviation of
% performance. The function also optionally returns the expected
% information gain ig from the comparison [a,b].
function [suggestion,ig]=optimal_pairwise(post_mean,post_cov,sigma)
    IGs=all_I_ij(post_mean,post_cov,sigma);
    [suggestion,ig]=argmax_mat(IGs);
end
    
    