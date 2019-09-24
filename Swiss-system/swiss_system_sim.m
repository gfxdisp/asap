function [q_sim,C_swiss] = swiss_system_sim(q_mat,number_of_observers,number_of_total_rounds, number_rand_rounds,N)
    % Condition of not comparing the same images more than once is violated
    number_of_swiss_rounds = number_of_total_rounds - number_rand_rounds;


    number_of_pairs = N/2;

    if (rem(N,2)==1)
        number_of_pairs = (N-1)/2;
        disp('Odd number of conditions, leaving');
    end

    C_swiss = zeros(N,N);

    % Go through the number of observers 
    for kk = 1:number_of_observers

        % Create matrix of the results for an observer
        C = zeros(N,N);
        % Go through random rounds
        for ii = 1:number_rand_rounds

            % Pairs 
            pairs = zeros(number_of_pairs,2);

            % Conditions from which to choose
            conds = [1:N];

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

                    C = simulate_observer_choice(q_mat,pairs(jj,1),pairs(jj,2),C);
            end

        end


        % Go through swiss sampling
        for ii = 1:number_of_swiss_rounds

            % array of pairs
            pairs = zeros(number_of_pairs,2);

            % Have an array of scores 1st col is id of the conditions
            % 2nd col is the current score of an image
            scores = zeros(N,2);
            scores(:,1) = [1:N]';
            scores(:,2) = sum(C,2);

            % Sort ascending based on the score
            scores= sortrows(scores,2);

            % Go through the number of pairs
            for jj=1:number_of_pairs
                
                number_conds = length(scores);

                rr = randi([1 number_conds],1,1);
                cc = randi([1 number_conds],1,1);
                range = 0;

                % if chosen the first condition
                if (rr == 1)
                    range = scores((rr+1),2)-scores(rr,2);
                % if last condition
                elseif (rr == number_conds)
                    range = scores(rr,2) - scores((rr-1),2);
                else
                    range = min(scores((rr+1),2) - scores(rr,2),...
                                scores(rr,2) - scores((rr-1),2));
                end

                % could have chosen the neighbouring condition, however, want to make
                % sure that if there are conditions with the same scores one is chosen
                % at random
                while (rr==cc) || abs(scores(rr,2)-scores(cc,2))>range
                    cc = randi([1 number_conds],1,1);
                end
                
                % set pairs and empty the score array with elements
                pairs(jj,1) = scores(rr,1);
                pairs(jj,2) = scores(cc,1);

                empty_els = [rr,cc];
                scores(empty_els,:) = [];

                C = simulate_observer_choice(q_mat,pairs(jj,1),pairs(jj,2),C);
                
            end

        end
        C_swiss = C_swiss + C;
    end

    q_sim = (sum(C_swiss,2)/number_of_observers)';
end
