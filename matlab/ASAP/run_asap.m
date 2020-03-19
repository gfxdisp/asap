function [pairs_to_compare] = run_asap(M,mode)

    assert(size(M,1)==size(M,2),'Matrix M must be square.')
    assert(size(M,1)>2,'More than two conditions are required for active sampling.')
    assert(strcmp(mode, 'mst') || strcmp(mode,'seq'),'Please pass valid mode, either seq or mst.')
    
    
    N = size(M,1);    
    G = unroll_mat(M);
    
    init_data.n_iter = 4;
    init_data.Ms=zeros(N,1);
    init_data.Vs=0.5*ones(N,1);
    init_data.Mgs=zeros(1,2); 
    init_data.Pgs=zeros(1,2); 
    init_data.Nc = N;
    init_data.G = G;
    init_data.prob_cmps = ones(N);
    
    [inf_mat,pairs_to_compare]=compute_information_gain_mat(N,init_data);
    if strcmp(mode,'mst')
        pairs_to_compare = compute_minimum_spanning_tree(inf_mat);
    end

end