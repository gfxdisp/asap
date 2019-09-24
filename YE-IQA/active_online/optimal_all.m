% Computes the optimal action (where the action can be either a pairwise
% comparison or absolute judgement) for a given multivariate normal
% posterior with mean post_mean, covariance post_cov, and absolute
% judgement cutoffs gamma_hat. sigma denotes the standard deviation of
% performance, cost_ratio denotes how many times more expensive it is to 
% perform an absolute judgement than a pairwise comparison.
% i.e. Cost of absolute judgement=cost_ratio*Cost of pairwise comparison.
%
% abs_or_pair is 0 if the optimal action is an absolute judgement, 1 if
% it's a pairwise judgement. suggestion is a scalar if abs_or_pair==0, and
% a vector in the form [a,b] (where the optimal comparison is a with b)
% if abs_or_pair==1.
function [abs_or_pair, suggestion]=optimal_all(post_mean,post_cov,gamma_hat,sigma,cost_ratio)
    [sugg_abs,ig_abs]=optimal_absolute(post_mean,post_cov,gamma_hat,sigma);
    [sugg_pair,ig_pair]=optimal_pairwise(post_mean,post_cov,sigma);
    w_ig_abs=ig_abs/cost_ratio;
    if w_ig_abs>ig_pair
        abs_or_pair=0;
        suggestion=sugg_abs;
    else
        abs_or_pair=1;
        suggestion=sugg_pair;
    end
end