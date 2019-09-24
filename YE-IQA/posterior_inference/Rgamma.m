% Evaluates the function \mathcal{R}_\gamma(s) defined in the paper, i.e.
% the Hessian (with respect to s) of \mathcal{F}_\gamma(s).
%
% We're  currently using numerical differentiation to do this (the
% derivest library). We could also (with some work) use automatic
% differentiation (autodiff), which would probably be faster, but this
% is a bit more complicated to implement for our code.
function R=Rgamma(s,mu,omega,gamma,P,M,sigma)
    f = @(x) Fgamma(x,mu,omega,gamma,P,M,sigma);
    
    %R = quick_ahess(f,s);% Automatic differentiation.
    R = hessian(f,s);% Numerical differentiation.
    %R = num_hess(f,s,eps^(1/3));% Fast numerical differentiation (less accurate).
end
