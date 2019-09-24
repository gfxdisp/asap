function [score, try_result]=init_score(data, mu, sigma, alpha, beta, para) 
    
    n_data=size(data,1);
    try_result=cell(n_data,2);
    gamma=getOpt(para,'gamma',0);    
    
    if length(unique(mu))==1 && length(unique(sigma))==1 && length(unique(alpha))==1 && length(unique(beta))==1
        
        [mu1_new, mu2_new, sigma1_new, sigma2_new, alpha_new, beta_new, KL_o, KL_a, win_prob]=online_update(mu(1),mu(1), sigma(1), sigma(1), alpha(1), beta(1), para);
        score=(KL_o+gamma*KL_a).*ones(n_data,1);        
        tmp=struct('mu1', mu1_new, 'mu2', mu2_new', 'sigma1', sigma1_new, 'sigma2', sigma2_new, 'alpha', alpha_new, 'beta', beta_new);
        for r=1:n_data
            try_result{r,1}=tmp;
            try_result{r,2}=tmp;
        end
        
    else
        
        score=zeros(n_data, 1);
        for r=1:n_data
            i=data(r,2);
            j=data(r,3);
            k=data(r,1);
            [try_result{r,1}.mu1, try_result{r,1}.mu2, try_result{r,1}.sigma1, try_result{r,1}.simga2, try_result{r,1}.alpha,  try_result{r,1}.beta,...
                KL_win_o, KL_win_a, win_prob]=online_update(mu(i), mu(j), sigma(i), sigma(j), alpha(k), beta(k), para);
            [try_result{r,2}.mu1, try_result{r,2}.mu2, try_result{r,2}.sigma1, try_result{r,2}.simga2, try_result{r,2}.alpha,  try_result{r,2}.beta,...
                KL_lose_o, KL_lose_a, lose_prob]=online_update(mu(j), mu(i), sigma(j), sigma(i), alpha(k), beta(k), para);
            score(r)=win_prob*(KL_win_o+gamma*KL_win_a)+lose_prob*(KL_lose_o+gamma*KL_lose_a);        
        end
        
    end

end