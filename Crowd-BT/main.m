%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation drive file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic inputs
clear 

n_anno = 8;  % number of workers
n_obj = 10;   % number of items
budget = 30;      % total available budget 
trials = 3;  % number of independent trials to run
u = 4; % Beta(u,v) is the generating distribution of workers' reliability
v = 1;

alpha0 = 4; % Beta(alpha0,beta0) is the prior of workers' reliability
beta0 = 1;
gamma = 1; % exploration-exploitaion tradeoff

%% Test

record = zeros(trials,budget); % record of ranking accuracy

for tr = 1:trials
    
    theta = drchrnd(ones(1,n_obj),1); % generate items
    rho = betarnd(u,v,1,n_anno); % generate workers
    
    % create data
    all_comb = combnk(1:n_obj, 2);
    [m,~] = size(all_comb);
    data = zeros(m*n_anno,3);
    for w = 1:n_anno
        before = (w-1)*m+1;
        after = w*m;
        data(before:after,1) = w;
        data(before:after,2:3) = all_comb;
    end
    
    % initial parameters
    mu = zeros(n_obj,1);
    sigma = ones(n_obj,1);
    alpha = alpha0 .* ones(n_anno,1);
    beta = beta0 .* ones(n_anno,1);
    
    active_para = struct('c', 0.1, 'kappa', 1e-4, 'taylor', 0, 'gamma', gamma, 'calc_iter', 1, ...
                        'sel_method', 'greedy', 'anno_threshold', 1e-4);

    [mu, sigma, alpha, beta, accuracy, hist]...
        = active_learning(data, budget, mu, sigma, alpha, beta, theta, rho, tr, active_para);
    
    record(tr,:) = accuracy;

end

%% Save simulation result

averaged = mean(record);
save(['./active_', active_para.sel_method,'_',num2str(n_anno),'_workers_',num2str(n_obj),'_items.mat'],...
    'averaged');
