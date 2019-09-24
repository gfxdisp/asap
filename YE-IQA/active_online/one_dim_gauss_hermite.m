% Approximates \int_{-\infty}^\infty e^{-x^2}f(x)\,\mathrm{d}x using
% Gauss-Hermite quadrature.
function integ=one_dim_gauss_hermite(f)
    N=50;
    [X, W] = hermquad(N);
    integ=0.0;
    threshw=max(W)*1e-6;
    for i=1:N
        wi=W(i);
        if wi>threshw
            integ=integ+wi*f(X(i));
        end
    end
end