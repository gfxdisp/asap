function [mu1_new, mu2_new, sigma1_new, sigma2_new, alpha_new, beta_new,...
    KL_o, KL_a, win_prob]=online_update(mu1, mu2, sigma1, sigma2, alpha, beta, para)

    c=getOpt(para, 'c', 25/6);
    kappa=getOpt(para, 'kappa', 1e-4);
    taylor=getOpt(para, 'taylor', false);

    tau=sqrt(sigma1^2+sigma2^2+2*c^2);
    e1=exp(mu1/tau);
    e2=exp(mu2/tau);
    e12=e1+e2;
    ae1=alpha*e1;
    be2=beta*e2;
    abe12=ae1+be2;
    
    win_prob=(ae1+be2)/(e12*(alpha+beta));
    
    tmp1=ae1/abe12-e1/e12;
    mu1_new=mu1+sigma1^2/tau*tmp1;
    mu2_new=mu2-sigma2^2/tau*tmp1;

    tmp2=(ae1*be2/(abe12^2)-e1*e2/(e12^2));
    sigma1_new = max(1+sigma1^2/tau^2*tmp2, kappa);
    sigma2_new = max(1+sigma2^2/tau^2*tmp2, kappa);
    
    KL_o=log(sigma1/sigma1_new)+(sigma1_new^2+(mu1_new-mu1)^2)/(2*sigma1^2)-0.5+...
         log(sigma2/sigma2_new)+(sigma2_new^2+(mu2_new-mu2)^2)/(2*sigma2^2)-0.5;
    
    if taylor
        C1=min(e1/e12+(sigma1^2/tau^2+sigma2^2/tau^2)*(e1*e2*(e2-e1)/(e12^3)), 1-kappa);
    else
        C1=e1/e12;
    end
    C2=1-C1;
    
    C=(C1*alpha+C2*beta)/(alpha+beta);
    m=(C1*(alpha+1)+C2*beta)*alpha/(C*(alpha+beta+1)*(alpha+beta));
    v=(C1*(alpha+2)+C2*beta)*(alpha+1)*alpha/(C*(alpha+beta+2)*(alpha+beta+1)*(alpha+beta));
     
    alpha_new=(m-v)*m/(v-m^2);
    beta_new=(m-v)*(1-m)/(v-m^2);
    KL_a=log(beta_func(alpha,beta))-log(beta_func(alpha_new, beta_new))-(alpha-alpha_new)*psi(alpha_new)...
         -(beta-beta_new)*psi(beta_new)+(alpha-alpha_new+beta-beta_new)*psi(alpha_new+beta_new);

end

function y = beta_func(z,w)
if nargin<2,
  error(message('MATLAB:beta:NotEnoughInputs'));
end
y = exp(gammaln(z)+gammaln(w)-gammaln(z+w));
end
