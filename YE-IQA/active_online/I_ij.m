% This function evaluates I(\mathcal{E}_{ij},Pr(\theta)) for a given i,j.
% It does this by splitting the calculation into three parts: p1,p2,p3.
%    p1 is E_\theta(p_{ij}\log(p_{ij})+q_{ij}\log(q_{ij})).
%    p2 is E_\theta(p_{ij})\log(E_\theta(p_{ij})).
%    p3 is E_\theta(q_{ij})\log(E_\theta(q_{ij})).
function IG=I_ij(mu,omega,sigma)
    integrand1=@(theta1,theta2) Eij_big_integrand([theta1,theta2],sigma);
    p1=two_dim_Etheta(mu,omega,integrand1);
    
    pij=@(theta1,theta2) pij_given_theta([theta1,theta2],sigma);
    qij=@(theta1,theta2) qij_given_theta([theta1,theta2],sigma);
    E_theta_pij=two_dim_Etheta(mu,omega,pij);
    E_theta_qij=two_dim_Etheta(mu,omega,qij);
    p2=E_theta_pij*log(E_theta_pij);
    p3=E_theta_qij*log(E_theta_qij);
    
    IG=p1-p2-p3;
end