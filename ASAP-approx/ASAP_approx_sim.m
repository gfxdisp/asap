function [M] = ASAP_approx_sim(q,n_comp,mode,M)

    N = size(M,1);
    G = [];
    comp= sum(M(:));
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

    while comp<n_comp
        [inf_mat,init_data,arg]=KLD_mat_approx(N,init_data);

        if strcmp(mode,'seq')
            G = simulate_observer_choice(q,arg(1),arg(2),G);
            comp=comp+1;
        else
            inf_mat = inf_mat+inf_mat';
            inf_mat = 1./inf_mat;
            inf_mat(inf_mat<0) = Inf;
            GrMST = graph(inf_mat);

            t = minspantree(GrMST);

            t_edges = sortrows(t.Edges,2);
            comparisons_to_perform = t_edges.EndNodes;
            for ii = 1:(N-1)
                G = simulate_observer_choice(q,comparisons_to_perform(ii,1),comparisons_to_perform(ii,2),G);

                comp=comp+1;
                if comp>=n_comp
                    break
                end

            end
        
            init_data.Pgs = [init_data.Pgs; zeros(N-2,2)];
            init_data.Mgs = [init_data.Mgs; zeros(N-2,2)];    
        end
        init_data.G = G;
    end
    
	for ii = 1:size(G,1)
        M(G(ii,1),G(ii,2)) = M(G(ii,1),G(ii,2))+1;
    end
end
