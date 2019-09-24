function [b, stats,pfit] = paired_comparisons(pcm, link_fcn)
% Given the pcm matrix of paired comparisons the corresponding scale values
% are returned
%
% Input:
%  pcm - paired comparison matrix 
%       (row m column n counts the number of times stimulus m is preferred
%       or selected over stimulus n)
%  link_fcn - link function to use for generalized linear model (GLM)
%       'logit' - for Bradley-Terry Model
%       'probit' - for Thurstone-Mosteller Model
%
% Output:
%   b - scale values whose order conforms with the columns/rows of pcm
%   stats - statistics based on the GLM fit

if nargin<2
    link_fcn = 'probit';
    fprintf(1,'Using Default Thurstone-Mosteller Model\n');
end


[M,N] = size(pcm);
X = zeros(N*(N-1)/2, N);

weights = zeros(N*(N-1)/2,1);
counts = zeros(N*(N-1)/2,1);

% Generate the model matrix X and the response vector Y (weights and
% counts)
k = 1;
for m = 1:(M-1)
    for n = (m+1):M
        weights(k) = pcm(m,n);
        counts(k) = pcm(m,n)+pcm(n,m);
        X(k,m) = 1;
        X(k,n) = -1;
        k = k + 1;
    end
end

%--------
% 

%--------
% Compute the Scale Values
% Let model have the form g(E{y}) = Xb, where 
%    link function is g(.)  
%    expectation operator is E{.}
%    X is the model matrix
%    b are the desired scale values


dist_y = 'binomial'; % distribution of the responses (y)
% Parameter Used for Bradley-Terry Model
% link_fcn = 'logit'; % link function

% Parameter Used for Thurstone-Mosteller Model
% link_fcn = 'probit';

% The constant term is set to zero (i.e., 'constant','off')
[b, dev, stats] = glmfit(X,[weights counts], dist_y ,'link', link_fcn,'constant','off');
%[b, dev, stats] = glmfit(X,[weights counts], dist_y, link_fcn);
%[b, dev, stats] = glmfit(X, [weights counts], dist_y);
% Offset = zeros(size(counts))
% Pwts = ones(size(counts))
% [b, dev, stats] = glmfit(X, [weights counts], dist_y, link_fcn, 'on', Offset, Pwts, 'off');
% [yhat, dylo, dyhi]=glmval(b, X, link_fcn, stats,'constant','off');
% ci=1.96.*stats.se./sqrt(counts(1));
% save('confidence_interval','ci');
%fprintf(1,'\tDeviance = %f\n', dev);
%fprintf(1,'\tDegrees of Freedom = %d\n', stats.dfe);


% Test the "fit" of the model to the data using the Deviance
pfit = 1-chi2cdf(dev, stats.dfe);
if pfit < 0.05
    fprintf(1,'\t***Fit is Suspect (p=%f)***\n',pfit);
    fprintf(1,'\tFor %d df, the deviance should be less than %f\n', stats.dfe, chi2inv(0.95, stats.dfe));
end

stats.dev = dev;
