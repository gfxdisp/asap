function kl = norm_kl(mu1,mu2,v1,v2)% KL divergence between two multivariate Gaussians (with diagonal covariance matrices).
    d = length(mu1);
    v2inv = 1 ./ v2;
    mdiffsq = (mu2-mu1).^2;
    total = sum(log(v2)) - sum(log(v1));
    total = total - d;
    total = total + sum(v1 .* v2inv);
    total = total + dot(v2inv, mdiffsq);
    kl = 0.5*total;
end