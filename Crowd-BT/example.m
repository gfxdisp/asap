
    
    number_conditions = 10;   
    number_comparisons = 100;      

    M = zeros(number_conditions);
    
    alpha0 = 5; % Beta(alpha0,beta0) is the prior of workers' reliability
    beta0 = 1;
    
    % A larger gamma will give more weight to the KL divergence terms 
    % related to annotator quality in the objective in Eq. (10), which 
    % means that we are willing to spend more to explore the quality of 
    % annotators. On the other hand, a smaller γ will result in relatively
    % more emphasis on exploiting the information in the observed pairwise 
    % comparisons. When there are gold samples, one could choose γ by 
    % testing the performance on the gold samples. However, when gold 
    % samples are not available, according to our experience, any 
    % γ ∈ [5, 10] could lead to much better performance than setting γ = 0 
    % (i.e., traditional active learning without exploring annotator 
    % quality) or γ = 1 (i.e., traditional information gain defined by KL
    % divergence)
    
    gamma = 0; % exploration-exploitaion tradeoff
    
    % create data
    all_comb = combnk(1:number_conditions, 2);
    [m,~] = size(all_comb);
    data = zeros(m,3);
    for w = 1:1
        before = (w-1)*m+1;
        after = w*m;
        data(before:after,1) = w;
        data(before:after,2:3) = all_comb;
    end

    % initial parameters
    mu = zeros(number_conditions,1);
    sigma = ones(number_conditions,1);
    alpha = alpha0 .* ones(1,1);
    beta = beta0 .* ones(1,1);

    active_para = struct('c', 0.1, ...
                         'kappa', 1e-4,...
                         'taylor', 0, ...
                         'gamma', gamma, ...
                         'calc_iter', 1, ...
                         'sel_method', 'multinomial', ...
                         'anno_threshold', 1e-4);

    [M] = active_learning(data, number_comparisons, mu, sigma, alpha, beta, theta, active_para,M);