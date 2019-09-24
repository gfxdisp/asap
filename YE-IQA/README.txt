This code is an online learning (non-batch) version of the following paper:
Ye, Peng and David Doermann. “Combining preference and absolute judgements in a crowd-sourced setting.” (2013). 
Terms, variables, and features mentioned in the rest of this README are discussed in more detail in the paper. Note that add_subfolders.m (in the main directory) must be run first in order to to add all the subfolders (containing important functions) to the MATLAB path.

The code has two main functions:

1. The function compute_posterior(mu,omega,M,sigma,P,tau) (in the posterior_inference folder) takes the following arguments:

mu, the mean of the (multivariate Gaussian) prior over the qualities q (denoted by s in the paper).
omega, the covariance of the (multivariate Gaussian) prior over q.
M, the absolute judgement matrix.
sigma, the standard deviation of condition performance (note that the paper and code implicitly assumes that the standard deviation of performance is the same for both pairwise comparisons and absolute judgements, although empirically they are not quite the same).
P, the pairwise comparison matrix.
tau, a small (a value of 0.01 or 0.005 is suggested) constant which is used to augment P and M in order to ensure unique minimizers for various functions.

It outputs an array [post_mean,post_cov,gamma_hat], where:

post_mean is the mean of the (multivariate Gaussian) posterior over the qualities q (denoted by s in the paper).
post_cov, the covariance of the (multivariate Gaussian) posterior over q.
gamma_hat, the maximum likelihood estimate of the absolute judgement cutoffs gamma.

The functions compute_posterior_absolute_only(mu,omega,M,sigma,tau) and compute_posterior_pairwise_only(mu,omega,sigma,P,tau) are identical to compute_posterior(mu,omega,M,sigma,P,tau) except that they take only matrices of absolute judgements and comparisons, respectively.



2. The function optimal_all(post_mean,post_cov,gamma_hat,sigma,cost_ratio) (in the active_online folder), which takes the following arguments:

post_mean, the mean of the (multivariate Gaussian) posterior over the qualities q (denoted by s in the paper).
post_cov, the covariance of the (multivariate Gaussian) posterior over q.
gamma_hat, the estimate for the absolute judgement cutoffs gamma.
sigma, the standard deviation of condition performance (note that the paper and code implicitly assumes that the standard deviation of performance is the same for both pairwise comparisons and absolute judgements, although empirically they are not quite the same).
cost_ratio, which denotes how many times more expensive it is to perform an absolute judgement than a pairwise comparison.
            i.e. Cost of absolute judgement=cost_ratio*Cost of pairwise comparison.

It outputs an array [abs_or_pair, suggestion], where:

abs_or_pair is 0 if the optimal action (under an expected information gain criterion) is an absolute judgement, 1 if it's a pairwise judgement. suggestion is a scalar if abs_or_pair==0, and a vector in the form [a,b] (where the optimal comparison is a with b) if abs_or_pair==1.

The functions optimal_absolute(post_mean,post_cov,gamma_hat,sigma) and optimal_pairwise(post_mean,post_cov,sigma) are similar, except that they return an array [suggestion,ig] (where suggestion is the optimal action with the same format as optimal_all, and ig is the expected information gain from the optimal action) for cases where only absolute judgements are available and only pairwise comparisons are available, respectively.



The examples folder contains two files which show how these functions can be used:

1. posterior_example.m simulates some data for q=0:4 and uses compute_posterior(mu,omega,M,sigma,P,tau) to estimate the posterior.
2. optimal_example.m takes the posterior from posterior_example.m and outputs the optimal next action.