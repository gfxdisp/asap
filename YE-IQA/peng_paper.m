function [P] = peng_paper(q,n_comp)
    next_comp=[1,2];
    N = size(q,2);
    M=zeros(4,N);% Absolute judgement matrix.
    P=zeros(N);% Comparison result matrix.
    sigma = 1.4865;
    mu=N*rand(1,N);
    omega=diag((1*ones(1,N)).^2);
    tau=5e-1;% tau should be fairly small
    for comp=1:n_comp
        P = simulate_observer_choice(q,next_comp(1),next_comp(2),P);
        
        [post_mean,post_cov] = compute_posterior(mu,omega,M,sigma,P,tau);

        IGs_ij=all_I_ij(post_mean,post_cov,sigma);
        next_comp=argmax_mat(IGs_ij);
    end
end