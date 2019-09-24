function [M] = HR_active_sim(q, n_cmps, M)


N= size(M,1);
totalnum=n_cmps;
lambda = 1e-1;
solve_model = 3;
% 1:uniform;  2:Bradley-Terry  3:Thurstone-Mosteller   4:Angular transform

[mapping_idx] = id_to_idxes(N);

i=[];
j=[];
for k=2:N
    j = [j;ones(k-1,1)*k];
    i = [i;(1:(k-1))'];
end

m = N*(N-1)/2;
d = sparse([1:m;1:m]',[i;j],[ones(1,m),-ones(1,m)],m,N);

w = zeros(m,1); 
pi = w;  


Sigma_0_inv=pinv(d'*((w*ones(1,N)).*d) + lambda*eye(N));
s_initial = Sigma_0_inv * (d' * (w .* (2*pi-1)));
s_initial = s_initial - mean(s_initial);


for num=sum(sum(M)):totalnum
    p = score2prob(s_initial, solve_model); 
    ii = (i-1)*N + i;
    jj = (j-1)*N + j;
    ij = (j-1)*N + i;
    d_S0_inv_d = Sigma_0_inv(ii) + Sigma_0_inv(jj) - 2*Sigma_0_inv(ij);  
    d_s0 = d*s_initial;
    d_S1_inv_d = d_S0_inv_d - d_S0_inv_d.^2./(1+d_S0_inv_d);

    ds_Sigma0_ds = ((1 + d_s0).^2-4*d_s0.*p).*d_S0_inv_d./(1+d_S0_inv_d).^2;
    informationgain = ds_Sigma0_ds - log(1- d_S1_inv_d) - d_S1_inv_d;

    [~, index]=sort(informationgain,'descend');
    maxgain_index=index(1);  

    Sigma0_inv_d = Sigma_0_inv*d(maxgain_index,:)';
    
    idx_1 = mapping_idx(maxgain_index,1);
    idx_2 = mapping_idx(maxgain_index,2);
    M_old = M;
    M = simulate_observer_choice(q,idx_1,idx_2,M);
    tmp = M - M_old;
    if tmp(idx_1,idx_2) == 1
        active_score = s_initial + Sigma0_inv_d*(1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
    else
        active_score = s_initial + Sigma0_inv_d*(-1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
    end
    Sigma_0_inv = Sigma_0_inv - Sigma0_inv_d*Sigma0_inv_d'/(1+d_S0_inv_d(maxgain_index));
    s_initial = active_score;

end

 