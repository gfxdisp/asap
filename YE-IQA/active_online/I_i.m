% This function evaluates I(\mathcal{E}_i,Pr(\theta)) for a given i.
% It does this by splitting the calculation into two parts: p1 and p2.
%    p1 is E_\theta(\sum_{k=1}^K p_{ik}\log(p_{ik}))).
%    p2 is \sum_{k=1}^K E_\theta(p_{ik})E_\theta(\log(p_{ik})).
function IG=I_i(mu0,sigma0,gamma_hat,sigma)
    p1_integrand=@(theta) entropy_pik_given_theta(theta,gamma_hat,sigma);
    p1=one_dim_Etheta(mu0,sigma0,p1_integrand);
    
    p2=0.0;
    K=numel(gamma_hat)-1;
    for k=1:K
        pik=@(theta) pik_given_theta(k,theta,gamma_hat,sigma);
        E_pik=one_dim_Etheta(mu0,sigma0,pik);
        p2=p2+E_pik*log(E_pik);
    end
    
    IG=p1-p2;
end