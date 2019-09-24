function [pairs_to_compare] = ts_sampling(M)

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
    
    pairs_to_compare = [];
    for ii = 2:N
            
       [~,pairs_to_compare_ii]=fill_trueskillmat(N,init_data,ii);
       
       pairs_to_compare = [pairs_to_compare; pairs_to_compare_ii];
       
    end

end