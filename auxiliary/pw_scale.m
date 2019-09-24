function [Q] = pw_scale( D, options )
% Scaling method for pairwise comparisons, also for non-balanced
% (incomplete) designs.
%
% [Q] = pw_scale( D, options )
%
% D - NxN matrix with positive integers. D(i,j) = k means that the
%     condition i was better than j in k number of trials.
% options - a cell array with the options. Currently recognized options:
%	   'prior' - type of the distance prior in the available options are:
%                'none': do not use prior, 'bounded': , 'gaussian': the
%                normalised sum of probabilities of observing a difference 
%                for all compared pairs of conditions. Set to 'none' by default. 
% Q - the JOD value for each method. The difference of 1 corresponds to
%     75% of answers selecting one condition over another.
%
% The condition with index 1 (the first row in D) has the score fixed at value
% 0. Always put "reference" condition at index 1. The measurement error is
% the smallest for the conditions that are the closest to the reference
% condition.
%
% The method scaled the data by solving for maximum-likelihood-estimator
% explaining the collected data. The principle is similar to Bradley-Terry
% model, but Gaussian and not logistic function is used to model relation
% between probabilities and distances in the scaled space. Similar
% scaling was proposed in:
%
% Silverstein, D., & Farrell, J. (2001). Efficient method for paired
% comparison. Journal of Electronic Imaging, 10, 394-398. doi:10.1117/1.1344187
%
% However, the method contains a number of extensions improving robustness,
% eliminating bias and allowing computing confidence intervals (with
% pw_scale_bootstrp function).
%
% Author: Rafal Mantiuk

% Revision history
% 2016-03-19 - Fixed graph connectivity patch; Replaced UA weights with
%              a conditional prior
% 2017-09-13 - Refined the prior and code simplification

% All elements must be non-negative integers 
if any(isinf(D(:))) || any(floor(D(:)) ~= D(:)) || any(D(:)<0)
    error( 'Matrix of comparisons contains invalid inputs');
end

if( ~exist( 'options', 'var' ) )
    options = {};
end

opt = struct();

% We don't the prior by default
opt.prior = 'none';
for kk=1:2:length(options)
    if( ~isfield( opt, options{kk} ) )
        error( 'Unknown option %s', options{kk} );
    end
    opt.(options{kk}) = options{kk+1};
end
 

if( size(D,1) ~= size(D,2) )
    error( 'The comparison matrix must be square' );
end

% The number of compared conditions
N = size( D, 1 ); 

options = optimset( 'Display', 'off', 'LargeScale', 'off' );

D_sum = D + D';
Dt = D';
nnz_d = (D_sum)>0;

% number of pairs compared at least once
comp_made = sum(nnz_d(:));


f = @(x)exp_prob(x);


% The methods tend to be more robust if starting scores are 0
Q = fminunc( f, zeros(N-1,1), options );

% Add missing leading 0-score for the first condition (not optimized)
Q = cat( 1, 0, Q );



    function P = exp_prob( q_trunc)
        q = cat( 1, 0, q_trunc ); % Add the condition with index 1, which is fixed to 0
                        
        sigma_cdf = 1.4826; % for this sigma normal cummulative distrib is 0.75 @ 1
        Dd = repmat( q, [1 N] ) - repmat( q', [N 1] ); % Compute the distances
        Pd = normcdf( Dd, 0, sigma_cdf ); % and probabilities  

        % Compute likelihoods
        prob = Pd(nnz_d);
        p = prob.^D(nnz_d) .*(1-prob).^Dt(nnz_d);        
        
        % Compute prior
        switch opt.prior
            case 'gaussian'
                
                prior = zeros(comp_made,1);
                for zz=1:N
                    for hh=1:N

                        n = D_sum(zz,hh);
                        %If the comparison has been performed
                        if n>0
                            k = D_wUA(zz,hh);
                            % Compute the probability of each distance
                            % according to all our answers
                            aux = prob.^k .* (1-prob).^(n-k);
                            prior = prior + (aux/sum(aux));

                        end
                    end
                end
                % The mean likelihood per answer is our prior (i.e., we compute
                % the probability of observing a certain distance according to
                % the rest of the answers in our comparison matrix)
            
            case 'bounded'
                q_range = max(q)-min(q);
                n_e = q_range+1;
                prior = max( NUA(nnz_d), 1/n_e - abs(D(nnz_d))/n_e.^2 );
                
            case 'none'
                prior = ones(comp_made,1);
            otherwise
                error( 'Unknown prior option %s', opt.prior );
        end
        
        P = -sum( log( max( p, 1e-400)  ) + log( max( prior + 0.1, 1e-400) )  );

    end

    function node_gr = connected_comp( G, node_gr, node, group )
        if( node_gr(node) ~= 0 )
            return;
        end
        node_gr(node) = group;
        for nn=find(G(node,:))
            node_gr = connected_comp( G, node_gr, nn, group );
        end
    end

end