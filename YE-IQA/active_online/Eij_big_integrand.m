% Evaluates p_{ij}\log(p_{ij})+q_{ij}\log(q_{ij}), given \theta, where
% p_{ij}=Pr(x_{ij}=1|theta), q_{ij}=Pr(x_{ij}=0|theta).
function integ=Eij_big_integrand(theta,sigma)
    pij=pij_given_theta(theta,sigma);
    if (pij>1e-8) && (pij<1-1e-8)
        qij=1-pij;
        integ=pij*log(pij)+qij*log(qij);
    else
        integ=0;
    end
end