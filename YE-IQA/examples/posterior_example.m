% This file can be used for testing the posterior inference code. It
% simulates pairwise comparison and absolute judgement data for
% q=[0:(n-1])], and then uses the posterior inference code to derive an
% estimate for q (which can then be compared with its true value).

tic
n=5;
q=2*randn(1,n);
%q=0:(n-1);% True condition qualities.
sigma=1/norminv(0.75);% Standard deviation in condition performance for pairwise comparisons. 

% Generate the pairwise comparison matrix. This is the matrix P from the
% paper.
num_pairwise=200;% Number of pairwise comparison experiments.
P=zeros(n);% Comparison matrix. This is the matrix P from the paper.
for foo=1:num_pairwise
    comp=randperm(n,2);
    xa=randn(1)*sigma+q(comp(1));
    xb=randn(1)*sigma+q(comp(2));
    if xb>xa
        P(comp(2),comp(1))=P(comp(2),comp(1))+1;
    else
        P(comp(1),comp(2))=P(comp(1),comp(2))+1;
    end
end

% Generate the rating observation matrix M. This is the matrix M from the
% paper.
num_absolute=100;% Number of absolute judgements.
K=6;% Number of cutoffs.
true_gamma=[-1e6 0 n*sort(rand(1,K-2)) 1e6];% Cutoff values (\gamma in the paper).
M=zeros(K,n);% Judgement matrix. This is the matrix M from the paper.

% c is a constant which represents how much more variable absolute
% judgement ratings are than pairwise comparisons.
% Absolute judgement standard deviation = c*Pairwise comparison standard
% deviation. Note that the method introduced in the paper actually
% implicitly assumes that c==1, but we simulate with c\neq 1 to give a more
% accurate performance assessment.
c=1.3;% This value is based off of empirical data.
csigma=sigma*c;

for foo=1:num_absolute
    i=randi(n);
    x_r=randn(1)*csigma+q(i);
    k=discretize(x_r,true_gamma);
    M(k,i)=M(k,i)+1;
end

% Generate prior.
mu=n*rand(1,n);
omega=diag((1*ones(1,n)).^2);
tau=5e-3;% tau should be fairly small
toc

tic
[post_mean,post_cov,gamma_hat] = compute_posterior(mu,omega,M,sigma,P,tau);
toc


% Normalizing post_mean to a 0:(n-1) scale.
q_est=zscore(post_mean);
q_est=q_est-q_est(1);
q_est=q_est/mean(q_est(2:end)./(1:(length(q_est)-1)));

% Display the true value of q and its estimate for comparison.
q
q_est