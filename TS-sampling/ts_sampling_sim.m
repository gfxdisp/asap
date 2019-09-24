function [M] = ts_sampling_sim(q,n_comp, M)

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
        for ii = 2:N
            
            [init_data,arg]=fill_trueskillmat(N,init_data,ii);
            comp=comp+1;
            G = simulate_observer_choice(q,arg(1),arg(2),G);
            init_data.G = G;
            if comp>=n_comp
                break 
            end
        end

    end
    
    for ii = 1:size(G,1)
        M(G(ii,1),G(ii,2)) = M(G(ii,1),G(ii,2))+1;  
    end

end