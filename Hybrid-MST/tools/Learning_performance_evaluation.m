function res = Learning_performance_evaluation(pcm, gt_score, model)
% according to current pcm, calculate the estimated score and std
% calculating the performance of the estimated value with the ground truth

% model could be 'logit' -- bradley-terry
% or 'probit' -- thurstone
% Score size: 1xn
% std size: 1xn


if strcmp(model,'logit') || strcmp(model,'probit')
[b,stats,pfit] = paired_comparisons(pcm,model);
elseif strcmp(model,'hodgerank')
    %% the input of hodgerank is a two-column matrix, with each observation result
    %% where the first column represents the stimulus is chosen than the other
    b= Hodgerank(pcm);
    stats = {};
    pfit = {};
else
    error('no matched pc model');
end

% [CC, Kendall, Spearman]=Correlation_performance(b,Score');
% 
% 
% res.b = b;
% res.stats = stats;
% res.pfit = pfit;
% res.performance = [CC, Kendall, Spearman];

res.kendall=corr(b,gt_score,'type','Kendall');
res.plcc=corr(b,gt_score);
%[res.CC, res.MAE, res.RMSE, res.ROCC]=statistic_analysis(b,gt_score);
            
