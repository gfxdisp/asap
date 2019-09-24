% Computes Pr(x_i=k|\theta).
% The theta argument for this function is s_i.
function pik=pik_given_theta(k,theta,gamma_hat,sigma)
    pik=normcdf((gamma_hat(k+1)-theta)/sigma)-normcdf((gamma_hat(k)-theta)/sigma);
end