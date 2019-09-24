function [Pearson, Kendall, Spearman] = Correlation_performance(x,y)

if size(x) ~= size(y)
    error('the size should be the same');
end

Pearson = corr(x,y,'type','Pearson');
Kendall = corr(x,y,'type','Kendall');
Spearman = corr(x,y,'type','Spearman');
