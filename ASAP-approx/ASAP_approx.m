function [pairs_to_compare] = ASAP_approx(M, mode)

    N = size(M,1);
    G = [];
    for ii=1:N
        for jj=1:N
            G=[G; repmat([ii jj], M(ii,jj),1)];
            M(ii,jj)=0;
        end
    end
    
    init_data.n_iter = 4;
    init_data.Ms=zeros(N,1);
    init_data.Vs=0.5*ones(N,1);
    init_data.Mgs=zeros(1,2); 
    init_data.Pgs=zeros(1,2); 
    init_data.Nc = N;
    init_data.G = G;
    init_data.prob_cmps = ones(N);

    [inf_mat,pairs_to_compare]=KLD_mat_approx(N,init_data);
    if strcmp(mode,'mst')
        inf_mat = inf_mat+inf_mat';
        inf_mat = 1./inf_mat;
        inf_mat(inf_mat<0) = Inf;
        GrMST = graph(inf_mat);

        t = minspantree(GrMST);

        t_edges = sortrows(t.Edges,2);
        pairs_to_compare = t_edges.EndNodes;
    end

end