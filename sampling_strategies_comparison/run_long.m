addpath('./data');
addpath('../ASAP-approx');
addpath('../ASAP');
addpath('../auxiliary');
addpath('../Crowd-BT');
addpath('../HR-active');
addpath(genpath('../Hybrid-MST'));
load('data/D_peng.mat')

N = 20;
N_exps = 4;
scaling_proc = 'trueskill';
dda = [11,12];
range_v = [0,5];

%MAT = D_peng(1:N,1:N);

q_mat.q = [];
q_mat.MAT = [];%MAT;%

t = rand(N_exps,N)*(range_v(2)-range_v(1))+range_v(1);
t = t - repmat(min(t, [], 2),[1,N]);
t = sort(t,2);
q_mat.rand_range = t;


start = N*(N-1)/2;
step = N*(N-1)/2;
finish = N*(N-1)/2*10;

cmps_arr = round([start:step/2:finish]);

file_name = strcat('./results_rand_s/data/',...
                   scaling_proc,'/xxxxxx',...
                   num2str(N),'_',...
                   num2str(range_v(2)),'_',...
                   strrep( strrep(num2str(dda),'  ',' '),' ','_'));


%poolobj = parpool(4);

simulate_long(q_mat,cmps_arr,N_exps, dda,file_name, scaling_proc)

delete(poolobj)



