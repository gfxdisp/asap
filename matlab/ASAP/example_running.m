% Create matrix with pairwise comparisons
M = [0,1,2,3,1;
     1,0,2,3,1;
     1,2,0,3,1;
     1,2,3,0,1;
     1,2,3,1,0];
              
% Run active sampling
pairs_to_compare = run_asap(M, 'mst');

% Calling print
display(pairs_to_compare)