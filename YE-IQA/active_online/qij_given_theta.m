% Computes Pr(x_{ij}=0|theta).
% The theta argument for this function is [s_i,s_j].
function qij=qij_given_theta(theta,sigma)
    qij=1-pij_given_theta(theta,sigma);
end