function pcm = Initial_learning(pcm)
%% this is used for the initialization of the pair comparison test


% we adopt the initial method of Ye Peng, 'active sampling for subjective IQA'

pcm = 0.5.*ones(size(pcm)) - diag(diag(0.5.*ones(size(pcm))));
