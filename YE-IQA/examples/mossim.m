niters=40;
n=5;% Number of conditions.
K=8;% Number of cutoffs.
sigma=1/norminv(0.75);% Standard deviation in condition performance for pairwise comparisons. 
cost_ratio=1;
c=1;% Empirical.
csigma=sigma*c;

trials=2;
total_moss=zeros(1,trials);
total_pcs=zeros(1,trials);
for trial=1:trials
    q=0:(n-1);% True condition qualities.
    true_gamma=[-1e6 0 4.5*sort(rand(1,K-2)) 1e6];% Cutoff values (\gamma in the paper).
    M=zeros(K,n);% Judgement matrix. This is the matrix M from the paper.
    P=zeros(n);% Comparison matrix. This is the matrix P from the paper.
    
    % Generate prior.
    mu=5*rand(1,n);
    omega=diag((2*ones(1,n)).^2);
    gamma_hat=[-1e6 linspace(0,(n+0.385),K-1) 1e6];
    tau=5e-3;% tau should be fairly small
    total_mos=0;total_pc=0;
    for iter=1:niters
        [abs_or_pair, suggestion]=optimal_all(mu,omega,gamma_hat,sigma,cost_ratio);
        if abs_or_pair==0% Absolute judgement.
            qi=q(suggestion);
            qi_samp=qi+randn()*csigma;
            judg=discretize(qi_samp,true_gamma);
            M(judg,suggestion)=M(judg,suggestion)+1;
            total_mos=total_mos+1;
            [iter, 0]
        else
            sugg1=suggestion(1);sugg2=suggestion(2);
            q1=q(sugg1);q2=q(sugg2);
            prob = normcdf(q1-q2, 0, sigma);
            if rand() < prob
                 P(sugg1,sugg2)=P(sugg1,sugg2)+1;
            else
                 P(sugg2,sugg1)=P(sugg2,sugg1)+1;
            end
            total_pc=total_pc+1;
            [iter,1]
        end
        [mu,omega,gamma_hat] = compute_posterior(mu,omega,M,sigma,P,tau);
    end
    total_moss(trial)=total_mos;
    total_pcs(trial)=total_pc;
end
total_moss
total_pcs

