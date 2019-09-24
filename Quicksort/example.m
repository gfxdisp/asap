number_comparisons = 100;
number_conditions = 10;

% Define and initialise quickM - a global variable with pairwise comparison
% matrix.
global quickM
quickM = zeros(number_conditions,number_conditions);

% Array with condition ranks
conditions_order=randperm(number_conditions);

% Run over allowed number of comparisons
while sum(quickM(:))<number_comparisons
    
    % Modify code in quicksort to add your experimental procedure to select
    % the outcome of the comparison
    conditions_order=quicksort(conditions_order,number_comparisons);
end

M = quickM;
