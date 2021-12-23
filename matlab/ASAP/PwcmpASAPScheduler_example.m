% Example of using PwcmpASAPScheduler

% We have 3 scenes (images or contents), each processed by two methods and
% each method is run using three different parameter values.
scene = { 'A', 'B', 'C' };
method = { 'methodA', 'methodB' };
param = [ 0, 1, 2 ];

% Create a table with all combinations of the elements from the sets
condition_table = create_factorial_table( scene, method, param );

% We want to compare relative quality within each scene: comapre all
% methods at all parameter values with each other, but without any
% cross-scene comparisons. 
sch = PwcmpASAPScheduler( 'test_results.csv', 'test_obs', condition_table, { 'scene' } );

% The code below simulates an experiment
N_batch = sch.get_pair_left(); % Get the number of conditions left in the current batch

for kk=1:N_batch % for each pairwise comparison in the batch
    
    [sch, stim_A, stim_B] = sch.get_next_pair(); % Get the next pair to compare
    
    fprintf( 1, 'Compare %d with %d\n', stim_A, stim_B );
    display( condition_table([stim_A stim_B],:) );
    
    % Simulate the response. In actual experiment, this is the place in which
    % the observer is asked to choose one of the two conditions.
    i_a = find( strcmp( condition_table.method{stim_A}, method ) );
    i_b = find( strcmp( condition_table.method{stim_B}, method ) );    
    p_sel = normcdf( i_b-i_a, 0, 2 );
    is_A_selected = rand()<=p_sel;
    
    % Record the result of the comparison. The result will be immediately
    % written to the file. 
    sch = sch.set_pair_result( is_A_selected );    
    
end
