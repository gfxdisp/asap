function p = score2prob(s,model)
N = length(s);
delta_s = s*ones(1,N)-ones(N,1)*s';
switch(model)
    case 1
        P = (delta_s + 1)/2;
    case 2
        P = 1./(1+exp(-delta_s));
    case 3
        P = (erf(delta_s) + 1)/2;
    case 4
        P = (sin(delta_s) + 1)/2;
    case 5 
        P = normcdf(delta_s,0,1.4826);
end
P = min(P,1);
P = max(P,0);
p = P(tril(ones(N))==0);