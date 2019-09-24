function [pairs] = swiss_system(C_init, mode)
% Function for swiss system sampling for pairwise comparisons.
% In the first rounds pairs of conditions are formed at random and later,
% conditions with similar scores compete together. 
% [pairs] = swiss_system(C_init, number_conditions, mode)
% C_init - matrix with pairwise comparisons collected so far (contains 0's
% if no comparisons have been collected)
% mode - 'rand' for random rounds and 'swiss' for rounds when similar
% condititions compete in pairs.
% The functions returns pairs - ids of elements to be compared in the
% experiment. If there is odd number of conditions, one (random) condition
% won't be compared.

    number_conditions = size(C_init,1);

    number_of_pairs = floor(number_conditions/2);
    
    % Go through random round
    if strcmp(mode, 'rand')

        % Pairs
        pairs = zeros(number_of_pairs,2);

        % Conditions from which to choose
        conds = [1:number_conditions];

        % Go through the number of pairs (N/2)
        for jj=1:number_of_pairs
            rr = 0;
            cc = 0;

            % Every time number of conditions to choose from is reduced
            new_N = length(conds);

            while rr==cc
                rr = randi([1 new_N],1,1);
                cc = randi([1 new_N],1,1);
            end

            % Set pairs
            pairs(jj,1) = conds(rr);
            pairs(jj,2) = conds(cc);

            % Remove chosen conditions
            empty_els = [rr,cc];
            conds(empty_els) = [];

        end
    else
        % array of pairs
        pairs = zeros(number_of_pairs,2);

        % Have an array of scores 1st col is id of the conditions
        % 2nd col is the current score of an image
        conds = zeros(number_conditions,2);
        conds(:,1) = [1:number_conditions]';
        conds(:,2) = sum(C_init,2);

        % Sort ascending based on the score
        conds= sortrows(conds,2);

        % Go through the number of pairs
        for jj=1:number_of_pairs

            number_conds = length(conds);

            rr = randi([1 number_conds],1,1);
            cc = randi([1 number_conds],1,1);
            range = 0;

            % if chosen the first condition
            if (rr == 1)
                range = conds((rr+1),2)-conds(rr,2);
            % if last condition
            elseif (rr == number_conds)
                range = conds(rr,2) - conds((rr-1),2);
            else
                range = min(conds((rr+1),2) - conds(rr,2),...
                    conds(rr,2) - conds((rr-1),2));
            end

            % could have chosen the neighbouring condition, however, want to make
            % sure that if there are conditions with the same scores one is chosen
            % at random
            while (rr==cc) || abs(conds(rr,2)-conds(cc,2))>range
                cc = randi([1 number_conds],1,1);
            end

            % set pairs and empty the score array with elements
            pairs(jj,1) = conds(rr,1);
            pairs(jj,2) = conds(cc,1);

            empty_els = [rr,cc];
            conds(empty_els,:) = [];

        end
    end
end

