% Evaluates \sum_{k=1}^K p_{ik}\log(p_{ik}) for a given \theta, where
% p_{ik}=Pr(x_i=k|\theta).
function entropy=entropy_pik_given_theta(theta,gamma_hat,sigma)
    K=numel(gamma_hat)-1;
    entropy=0.0;
    for k=1:K
        pik=pik_given_theta(k,theta,gamma_hat,sigma);
        % Need this if/else to ensure that entropy+=0 when pik=0.
        if pik>1e-8
            entropy=entropy+pik*log(pik);
        end
    end
end