% Evaluates E_\theta(f(\theta)), where \theta is bivariate normal with
% mean mu and covariance matrix omega. It does this by Gauss-Hermite
% quadrature. "A note on multivariate Gauss-Hermite quadrature" by
% Peter Jackel contains more information on this, and the R MultiGHQuad
% package has a more general implementation (in R).
function integ=two_dim_Etheta(mu,omega,f)
    N=50;
    [X, W]=hermquad(N);
    W=W/sqrt(pi);
    X=X*sqrt(2);

    Wmat=W*W.';
    [X1,X2] = meshgrid(X,X);
    Xmat = [X1(:) X2(:)].';

    
    [eigvects,diag_eigvals]=eig(omega);
    rotation=eigvects*sqrt(diag_eigvals);
    Xmat=rotation*Xmat+mu.';

    integ=0.0;
    Nsq=N^2;
    threshw=max(Wmat)*1e-6;
    for i=1:Nsq
        wi=Wmat(i);
        if wi>threshw
            integ=integ+wi*f(Xmat(1,i),Xmat(2,i));
        end
    end
end