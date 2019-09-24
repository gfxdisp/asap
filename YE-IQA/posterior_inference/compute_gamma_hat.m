% Computes \hat{\gamma}, the MLE for \gamma, which is the maximizer of
% log(Pr(M,P|\gamma)). The optimization iself is over gamma_part, which we
% define to be \gamma_2,\gamma_3,\cdots,\gamma_{\textsc{k}-1}.
% \gamma_0,\gamma_1, and \gamma_\textsc{k} are always set to be -\infty, 0,
% \infty, respectively.
function gamma_hat=compute_gamma_hat(P,M,sigma,mu,omega)
    K=size(M,1);
    n=size(M,2);
    % This default guess for gamma is used when either there's no absolute
    % judgement data, or the optimizer for computing \hat{\gamma} fails.
    default_gamma=[-1e6 0 linspace(0.5,(n-0.5),K-2) 1e6];
    if max(M(:))<1
        gamma_hat=default_gamma;
    else
        fun=@(gamma_part) -loglikeMP_part(P,M,gamma_part,sigma,mu,omega);% Objective function.

        % Define the linear constraints in the form A*gamma_part\leq b.
        diags=[ones(K-2,1) -ones(K-2,1)];
        A = full(spdiags(diags,[-1 0], K-1, K-2));
        b = [zeros(K-2,1); 1e6];
        % MATLAB fmincon requires equality constraints, so we just leave them
        % empty (since there are none).
        Aeq=[];beq=[];

        % We put an upper bound on the values of gamma, since otherwise the
        % optimizer starts evaluating very large values for some of the entries
        % in \gamma, which crashes the optimizer.
        lb=zeros(1,K-2);ub=n*ones(1,K-2);

        gamma_part0=sort(rand(1,K-2));% Initial guess.
        opts1=  optimset('display','off');% Suppress the MATLAB optimization console output.
        

        % -1e6 and 1e6 are just approximations for \pm\infty.
        try
            gamma_hat_part = fmincon(fun,gamma_part0,A,b,Aeq,beq,lb,ub,[],opts1);
            gamma_hat=[-1e6 0 gamma_hat_part 1e6];
        catch% Use default gamma if optimization fails.
            gamma_hat=default_gamma;
        end
    end
end