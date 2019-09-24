% Evaluates E_\theta(f(\theta)), where \theta is a scalar and is normally
% distributed with mean mu0 and covariance matrix sigma0. It does this by
% converting the expectation into an integral of the form
% \int_{-\infty}^\infty e^{-x^2}g(x)\,\mathrm{d}x, which can be
% approximated efficiently by Gauss-Hermite quadrature.
function integ=one_dim_Etheta(mu0,sigma0,f)
    sigma0sqrt2=sigma0*sqrt(2);
    f2=@(theta) f(sigma0sqrt2*theta+mu0);
    
    integ=one_dim_gauss_hermite(f2);
    integ=integ*pi^(-0.5);
end