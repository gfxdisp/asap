% Create matrix with pairwise comparisons
M = [0,1,2,3,1;
     1,0,2,3,1;
     1,2,0,3,1;
     1,2,3,0,1;
     1,2,3,1,0];
              
% Run active sampling
[pairs_to_compare, Mean, Std] = run_asap(M, 'mst');

% Calling print
display(pairs_to_compare)

% Converting to JODs
JOD = (Mean-mean(Mean))/std(Mean)*1.4826;
JOD = JOD-JOD(1);
JODstd = Std/std(Mean)*1.4826;