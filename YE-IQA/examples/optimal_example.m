% This file can be used to test the active_online code. It takes estimates
% for post_mean, post_cov, and gamma_hat from posterior_example, and
% outputs the optimal comparison/judgement for a given cost_ratio.
cost_ratio=1;
tic
[abs_or_pair, suggestion]=optimal_all(post_mean,post_cov,gamma_hat,sigma,cost_ratio)
toc