% Evaluates the prior log(Pr(s)). Note that we make this separate from the 
% logmvnpdf function to make it easier to add stuff to this function later
% if necessary.
function loglike=logprior(s,mu,omega)
    loglike=logmvnpdf(s,mu,omega);
end