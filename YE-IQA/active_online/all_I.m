% This function evaluates I(\mathcal{E}_i,Pr(\theta)) for all i and
% I(\mathcal{E}_{ij},Pr(\theta)) for all i,j. The results are returned as
% two separate arrays/matrices.
function [IGs_i, IGs_ij]=all_I(post_mean,post_cov,gamma,sigma)
    IGs_i=all_I_i(post_mean,post_cov,gamma,sigma);
    IGs_ij=all_I_ij(post_mean,post_cov,sigma);
end