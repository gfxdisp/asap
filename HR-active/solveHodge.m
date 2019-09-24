function [s,Sigma0] = solveHodge(d,w,pi,model,lambda)
N = size(d,2);
m = length(pi);
Sigma0 = d'*diag(w)*d + lambda*eye(N);
switch(model)
    case 1
        y = 2*pi - 1;
    case 2
        pi = max(pi,eps*ones(1,m));
        pi = min(pi,1-eps*ones(1,m));
        y = log(pi./(1-pi));
    case 3
        pi = max(pi,eps*ones(1,m));
        pi = min(pi,1-eps*ones(1,m));
        y = erfinv(2*pi-1);
    case 4
        y = asin(2*pi-1);
end
[s,flag]= lsqr(Sigma0,d'*diag(w)*y');
s = s - mean(s);



