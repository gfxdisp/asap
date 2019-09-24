function [mean_corr,mean_rmse]=simulate_exp_bootstrp_dif_designs(q_mat,design,n_comp,N_exps,scaling_proc)

    q_mat.q = q_mat.rand_range(1,:);
    N = max(length(q_mat.q),size(q_mat.MAT,1));
    q_s = zeros(N_exps,N);
    M_cmp = zeros(N,N);
    M = zeros(N,N);
    
    if size(q_mat.MAT,1)>0 
        q_true = pw_scale(q_mat.MAT)';
    end
    time_array = [];       
    for ee=1:N_exps
        if size(q_mat.MAT,1)==0
            q_mat_local = q_mat.rand_range(ee,:);
            q_true = q_mat.q;
        else
            q_mat_local = q_mat.MAT;
        end
        M = zeros(N,N);
        
        tic
        switch design
            case 1
                each = ceil(n_comp/(0.5*N*(N-1)));
                M = triu(ones(N,N),1).*each;
            case 2
                q_p = 1:N;
                M = zeros(N,N);
                each = ceil(n_comp/(N-1));
                for ii = 1:(N-1)
                    M(q_p(ii),q_p(ii+1)) = each;
                end
            case 3
                n_obs = ceil(n_comp/(4.5*N));
                M = swiss_pairing_sim(q_mat_local,n_obs,9,3, N);
            case 4
                [~,~,num_cmp_per_obs] = calc_numb_cmp(N);
                n_obs = round(n_comp/num_cmp_per_obs);
                n_rand_obs = round(0.3*n_obs);
                M = ard_sim(q_mat_local,n_obs,n_rand_obs,M);
            case 5
                M = peng_paper(q_mat_local,n_comp);
            case 6
                M = quicksort_sim(q_mat_local,N,n_comp);
            case 7
                M = ts_sampling_sim(q_mat_local,n_comp, M);
            case 8
                M = run_crowd_bt(q_mat_local,n_comp,M);
            case 9
                M = run_hodge_rank(q_mat_local, n_comp,M);
            case 10
                M = hybridMST(q_mat_local, n_comp, M);   
            case 11
                M = ASAP_sim(q_mat_local,n_comp,'mst',M);
            case 12
                M = ASAP_approx_sim(q_mat_local,n_comp,'mst',M);
        end
        res_t = toc;
        
        time_array = [time_array,res_t];
        M = remove_excess_and_simulate(design, M, q_mat_local, n_comp);
        
        M_cmp = M_cmp + ((M+M.')/2)/N_exps;
        if strcmp(scaling_proc,'trueskill') 
            [ms,vs]= ts_M(M);
        else
            ms = pw_scale(M);
            vs = zeros(1,N);
        end
        ms  =(ms-mean(ms))/std(ms);
        vs = vs/std(ms);
        q_s(ee,:) = ms;
        q_true  =(q_true-mean(q_true))/std(q_true);

        corr_arr(ee) = corr(q_s(ee,:)', q_true', 'Type', 'Spearman');
        rmse_arr(ee) = sqrt(mean((q_true - q_s(ee,:)).^2));
    end
    disp(['Time ms/cmp:', num2str(mean(time_array)/n_comp*1000)]);

    mean_corr = mean(corr_arr);
    mean_rmse = mean(rmse_arr);
end