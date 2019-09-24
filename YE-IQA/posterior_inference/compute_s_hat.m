% Computes the minimizer \hat{s}_\gamma of \mathcal{F}_\gamma(s) defined in
% the paper.
function s_hat=compute_s_hat(mu,omega,gamma,P,M,sigma)
    % s0 is the initial guess for the function minimizer.
    s0=rand(1,numel(mu));
    % Define anonymous function in order to pass extra parameters to the
    % objective function.
    f=@(s) Fgamma(s,mu,omega,gamma,P,M,sigma);
    opts1=  optimset('display','off');% Suppress the MATLAB optimization console output.
    s_hat=fminunc(f,s0,opts1);
end