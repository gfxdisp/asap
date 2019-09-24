total_rounds = 9;
random_rounds = 3;
ard_rounds = total_rounds - random_rounds;
number_conditions = 10;
C_init = zeros(number_conditions);

[~,~,numb_cmps] = calc_numb_cmp(number_conditions);
disp(['Total number of comparisons: ', num2str(numb_cmps*total_rounds)]);

% Go over random rounds
for ii = 1:random_rounds
    [pairs_to_compare] = ard(C_init, 'rand');
    for jj = 1:size(pairs_to_compare,1)
        pair = pairs_to_compare(jj,:);
        disp(['Rand round: ',num2str(ii),', observer compares: ',num2str(pair(1)),' ',num2str(pair(2))])
        
        % Simulate observer choice
        if rand()>0.5
            C_init(pair(1), pair(2)) = C_init(pair(1), pair(2))+1;
        else
            C_init(pair(2), pair(1)) = C_init(pair(2), pair(1))+1;
        end
    end
end

% Go over ard rounds
for ii = 1:ard_rounds
    [pairs_to_compare] = ard(C_init, 'ard');
    for jj = 1:size(pairs_to_compare,1)
        pair = pairs_to_compare(jj,:);
        disp(['ARD round: ',num2str(ii),', observer compares: ',num2str(pair(1)),' ',num2str(pair(2))])
        
        % Simulate observer choice
        if rand()>0.5
            C_init(pair(1), pair(2)) = C_init(pair(1), pair(2))+1;
        else
            C_init(pair(2), pair(1)) = C_init(pair(2), pair(1))+1;
        end
    end
end