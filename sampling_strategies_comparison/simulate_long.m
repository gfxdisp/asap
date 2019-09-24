function  simulate_long(q_mat,cmps_arr,N_exps, methods, file_name, scaling_proc)
    q_mat.q = q_mat.rand_range(1,:);
    
    N = max(numel(q_mat.q), size(q_mat.MAT,1));
    n_methods = numel(methods);
    var_names =  {'design', 'comp', 'ci', 'corr', 'rmse','mean_var', 'cmps_per_n_conds'};
    table_res = cell2table(cell(0,7), 'VariableNames',var_names);
    
    q_s = zeros(n_methods,size(cmps_arr,2),N_exps,N);
    v_s = zeros(n_methods,size(cmps_arr,2),N_exps,N);
    corr_res = zeros(1,N_exps);
    rmse_res = zeros(1,N_exps);
    var_res = zeros(1,N_exps);

    mm_c = 0;
    cc_c = 0;
    for mm = methods
        mm_c = mm_c+1;
        parfor ee=1:N_exps
            q_mat_local = q_mat.MAT;
            if size(q_mat.MAT,1)==0
                q_mat_local = q_mat.rand_range(ee,:);
            else
                q_mat_local = q_mat.MAT;
            end

            time_bootstrap = 0;
            M =zeros(N);
            cc_c = 0;
            tmp_ms_cc = zeros(size(cmps_arr,2), N);
            tmp_vs_cc = zeros(size(cmps_arr,2), N);
            for n_comp=cmps_arr
                cc_c = cc_c+1;
                disp(['mm: ',num2str(mm),' ee: ',num2str(ee),' cc: ',num2str(n_comp), ' st: ', num2str(n_comp/(N*(N-1)/2))])
                
                switch mm                    
                    case 7
                        [M] = trueskill_sampling(q_mat_local,n_comp, M);
                    case 8
                        M = run_crowd_bt(q_mat_local,n_comp,M);
                    case 9 
                        M = run_hodge_rank(q_mat_local, n_comp,M);
                    case 10
                        M = hybridMST(q_mat_local, n_comp, M);
                    case 11                
                        [M] = ASAP_sim(q_mat_local,n_comp,'mst',M);
                    case 12
                        [M] = ASAP_approx_sim(q_mat_local,n_comp,'mst',M);
                end
                %save(strcat('M_q_',num2str(cc_c)),'M','q_mat_local');

                if strcmp(scaling_proc,'trueskill') 
                    [ms,vs]= ts_M(M);
                else
                    ms = pw_scale(M);
                    vs = zeros(1,N);
                end
                ms  =(ms-mean(ms))/std(ms);
                vs = vs/std(ms);
                tmp_ms_cc(cc_c,:) =ms;
                tmp_vs_cc(cc_c,:) =vs;

            end
            
            q_s(mm_c,:,ee,:) = tmp_ms_cc;
            v_s(mm_c,:,ee,:) = tmp_vs_cc;
        end
        %save(strcat(file_name,'.mat'),'q_s','v_s','methods','q_mat')
        cc_c = 0;
        for n_comp = cmps_arr
            cc_c = cc_c+1;
            for ee=1:N_exps
            	if size(q_mat.MAT,1)==0
                    q_mat.rand_range(ee,:)= (q_mat.rand_range(ee,:)-mean(q_mat.rand_range(ee,:)))/std(q_mat.rand_range(ee,:));
                else
                    q_mat.rand_range(ee,:) = ts_M(q_mat.MAT);    
                    q_mat.rand_range(ee,:)= (q_mat.rand_range(ee,:)-mean(q_mat.rand_range(ee,:)))/std(q_mat.rand_range(ee,:));
                
                end
                corr_res(1,ee) = corr(squeeze(q_s(mm_c,cc_c,ee,:)), q_mat.rand_range(ee,:)', 'Type', 'Spearman');
                rmse_res(1,ee) = sqrt(mean((q_mat.rand_range(ee,:)' - squeeze(q_s(mm_c,cc_c,ee,:))).^2));
                var_res(1,ee) = mean(v_s(mm_c,cc_c,ee,:));
            end

            q_m = mean(squeeze(q_s(mm_c,cc_c,:,:)));
            e_up = prctile(q_s(mm_c,cc_c,:,:),97.5)-q_m;
            e_low = q_m - prctile(squeeze(q_s(mm_c,cc_c,:,:)),2.5);
            
            table_res = [table_res;cell2table(cell({mm, n_comp, mean(mean( cat( 2, e_up(2:end), e_low(2:end) ) )), mean(corr_res),mean(rmse_res),mean(var_res),n_comp/(N*(N-1)/2)}), 'VariableNames',var_names)];
            
        end
        writetable(table_res, strcat(file_name,'.csv'));
    end
end