% Computes Pr(x_{ij}=1|theta).
% The theta argument for this function is [s_i,s_j].
function pij=pij_given_theta(theta,sigma)
    pij=normcdf((theta(1)-theta(2))/(sigma*sqrt(2)));
end