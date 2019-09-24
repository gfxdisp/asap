addpath('./data');
addpath('../swiss_system');
addpath('../adaptive_rectangular_design');
addpath('../ASAP-approx');
addpath('../ASAP');
addpath('../auxiliary');
addpath('../Crowd-BT');
addpath('../HR-active');
addpath('../Quicksort');
addpath(genpath('../Hybrid-MST'));
addpath(genpath('./'));


N = 20;
N_exps = 10;
range_v = [0,5];
dda = [1,2,3,4,6];


q_mat.q = [];
q_mat.MAT = [];

t = rand(N_exps,N)*(range_v(2)-range_v(1))+range_v(1);
t = t - repmat(min(t, [], 2),[1,N]);
q_mat.rand_range = sort(t,2);

cmps_arr = [(N*(N-1)/2):(N*(N-1)/4):((N*(N-1)/2)*10)];


scaling_proc = 'trueskill';

file_name = strcat('./results_rand_s/data/',...
                   scaling_proc,'/tmpxxx',...
                   num2str(N),'_',...
                   num2str(range_v(2)),'_',...
                   strrep( strrep(num2str(dda),'  ',' '),' ','_'));
               
fh = fopen( file_name, 'w' );
fprintf( fh, 'design, comp, ci, corr, rmse, cmps_per_n_conds\n' );               

for dd=dda
    for cc=cmps_arr
        disp(['design type = ', num2str(dd),  'comparisons = ', num2str(cc)]);
        [mean_corr,mean_rmse] = simulate_exp_bootstrp_dif_designs(q_mat,dd,cc,N_exps,scaling_proc);
        fprintf( fh, '%d, %d, %g, %g, %g\n', dd, cc, mean_corr, mean_rmse, cc/(N*(N-1)/2));
     end
 end


fclose( fh );
