function [pairs_to_compare, Mean, Std] = run_asap(M,mode)
% Run active sampling for pairwise comaprisons method to find the next pair
% (or batch of pairs) to compare. 
%
% M - a square pairwise comparison matrix, where M(i,j)=k means that the
% condition i was selected as better than condition j k-times. 
% mode - 'mst' for the batch mode (minimum spanning tree) or 'seq' for the
%        sequential mode.

    assert(size(M,1)==size(M,2),'Matrix M must be square.')
    assert(size(M,1)>2,'More than two conditions are required for active sampling.')
    assert(strcmp(mode, 'mst') || strcmp(mode,'seq'),'Please pass valid mode, either seq or mst.')
            
    
    N = size(M,1);
    
    if sum(M(:))==0 && strcmp(mode, 'mst')
        % For the very first batch of comparisons, return a random set of
        % pairs
        rp = randperm(N);
        rord = randperm(N-1);
        pairs_to_compare = zeros(N-1,2);
        for kk=1:N-1
            ind = rord(kk);
            pairs_to_compare(ind,:) = rp(ind:(ind+1));
        end
        Mean = nan(N,1);
        Std = nan(N,1);
        return;
    end
    
    % turn matrix m into an array with dimensionality (numb_comparisons, 2), where in row kk 
    % tupple [3,4] means that in the kth comparison condition 3 was chosen over condition 4
    G = unroll_mat(M);
    
    % Set starting parameters for the TrueSkill
    init_data.n_iter = 4;
    init_data.Ms=zeros(N,1);
    init_data.Vs=0.5*ones(N,1);
    init_data.Mgs=zeros(1,2); 
    init_data.Pgs=zeros(1,2); 
    init_data.Nc = N;
    init_data.G = G;
    init_data.prob_cmps = ones(N);
    
    % Compute the information gain matrix and return a single pair maximizing information gain
    [inf_mat,init_data,pairs_to_compare]=compute_information_gain_mat(N,init_data);
    if strcmp(mode,'mst')
        pairs_to_compare = compute_minimum_spanning_tree(inf_mat);
    end
    
    Mean = init_data.Ms;
    Std = init_data.Vs;
    
end
