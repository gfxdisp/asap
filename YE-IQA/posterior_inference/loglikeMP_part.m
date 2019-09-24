% Define a version of the loglikeMP function (which evaluates
% log(Pr(M,P|\gamma))) that takes only \gamma_2,\gamma_3,\cdots,
% \gamma_{\textsc{k}-1} as input (and sets \gamma_0,\gamma_1, and
% \gamma_\textsc{k} are always set to be -\infty, 0, \infty, respectively.
% This function is needed to ensure that the function we're minimizing for
% compute_gamma_hat is in the correct format for the MATLAB optimizer.
function loglike=loglikeMP_part(P,M,gamma_part,sigma,mu,omega)
    % -1e6 and 1e6 are just approximations for \pm\infty.
    gamma=[-1e6 0 gamma_part 1e6];
    loglike=loglikeMP(P,M,gamma,sigma,mu,omega);
end