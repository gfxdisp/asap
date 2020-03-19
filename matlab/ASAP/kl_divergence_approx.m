function kl = kl_divergence_approx(mu1,mu2,v1,v2)
    total = sum(log(v2)) - sum(log(v1));
    total = total + sum(v1 ./ v2);
    total = total + dot(1 ./ v2, (mu2-mu1).^2);
    kl = total;
end