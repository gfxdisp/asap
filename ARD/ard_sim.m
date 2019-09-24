function [C] = ard_sim(q_mat,number_of_observers,number_rand_rounds,M)

    N = size(M,1);

    if (isprime(N))
        disp('Prime number of conditions, leaving');
        return
    end

    [numb_rows_vortex_M,numb_cols_vortex_M,~] = calc_numb_cmp(N);

    C = zeros(N,N);

    % Go through the number of observers 
    for oo = 1:number_of_observers

        M = zeros(numb_rows_vortex_M,numb_cols_vortex_M);
        
        % Have an array of scores 1st col is id of the conditions
        % 2nd col is the current score of an image
      	scores = zeros(N,2);
        scores(:,1) = [1:N]';
        
        % Go through random rounds
        if oo<=number_rand_rounds
            scores(:,1)= (randperm(length(scores(:,1))));
        else
            scores(:,2) = pw_scale(C,0);
            % Sort ascending based on the score
            scores = sortrows(scores,2);
        end
        
        [M] = recursive_filling(scores(:,1)',M,'right',1,1,numb_cols_vortex_M+1,numb_rows_vortex_M+1,0,0);
        pairs = construct_pairs(M);
        number_of_pairs = size(pairs,1);
        
        % Go through the number of pairs
        for jj=1:number_of_pairs
            C = simulate_observer_choice(q_mat,pairs(jj,1),pairs(jj,2),C);   
        end
        
    end

end


