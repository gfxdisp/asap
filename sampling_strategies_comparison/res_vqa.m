addpath('./data');
addpath('../swiss_system');
addpath('../active_ts_online_true_sc');
addpath('../active_ts_online_sc');
addpath('../auxiliary_functions');
addpath('../CrowdBT');
addpath('../hodge_rank');
addpath('../pwscale');
addpath('../quicksort');
addpath('../trueskill_sampling');
addpath(genpath('../Hybrid-MST'));

N = 16;
N_exps = 20;

q_mat.q = [];
range_v = [0,20];
t = rand(N_exps,N)*(range_v(2)-range_v(1))+range_v(1);
t = t - repmat(min(t, [], 2),[1,N]);
t = sort(t,2);
q_mat.rand_range = t;

start = N*(N-1)/2;
step = N*(N-1)/2;

dda = [31];
scaling_method = 'thurstone';
%poolobj = parpool(4);
for ff = 1:3
    load(strcat('data/VQA_Data_mat/', num2str(ff),'.mat'));
    q_mat.MAT = mat_to_save;
    full_cmps = sum(q_mat.MAT(:));
    cmps_arr = round([start:step/2:full_cmps/2]);
    
    res_file = strcat('./results_real_d/VQA/',scaling_method,'/',strrep( strrep(num2str(dda),'  ',' '),' ','_'),'_',num2str(ff));

    simulate_long(q_mat,cmps_arr,N_exps, dda,res_file, scaling_method)


end
delete(poolobj)
%% Types of design
% 1: full design
% 2: nearest conditions 
% 3: swiss design
% 4: adaptive squares
% 5: Peng Ye active sampling paper
% 6: Quicksort
%15: crowd bt
%16: hodge rank

% 7: TS - Pure Differential Entropy
% 8: TS - Pure KL-Divergence
% 9: TS - Differential Entropy/KL-Divergence Hybrid
%10: TS - Differential Entropy/KL-Divergence Hybrid + IG + weighted
%11: TS -  KL-Divergence MST
%12: TS -  Differential MST
%13: TS -  Differential+KL MST
%14: TS -  Differential+KL+IG MST
%17: TS - Weighted hybrid