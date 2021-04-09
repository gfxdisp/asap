function [pairs] = ard(C_init, mode)

    N = size(C_init,1);

    if (isprime(N))
        disp('Prime number of conditions, leaving');
        return
    end

    [numb_rows_vortex_M,numb_cols_vortex_M,~] = calc_numb_cmp(N);

    M = zeros(numb_rows_vortex_M,numb_cols_vortex_M);

    % Have an array of scores 1st col is id of the conditions
    % 2nd col is the current score of a condition
    scores = zeros(N,2);
    scores(:,1) = [1:N]';

    % Go through random rounds
    if strcmp(mode,'rand')
        scores(:,1)= (randperm(length(scores(:,1))));
    else
        scores(:,2) = sum(C_init,2);
        % Sort ascending based on the score
        scores = sortrows(scores,2);
    end

    [M] = recursive_filling(scores(:,1)',M,'right',1,1,numb_cols_vortex_M+1,numb_rows_vortex_M+1,0,0);
    pairs = construct_pairs(M);

end