addpath ('tools')

% optimisation procedure throws warnings, switch them off for speed
% reasons
warning('off','stats:glmfit:IllConditioned')
    
    
    
number_conditions  =  10;
number_comparisons = 100;
M = zeros(number_conditions);
    
    
%% test on mst model

trial = sum(sum(M));
trial_i = sum(sum(M));

full_pair_num = number_conditions*(number_conditions-1)/2;
ii = 0;

while ii < number_comparisons
    [pairs_to_compare] = hybrid_MST(M);

    for jj = 1:size(pairs_to_compare,1)
        pair = pairs_to_compare(jj,:);
        disp(['Round: ',num2str(ii),', observer compares: ',num2str(pair(1)),' ',num2str(pair(2))])
        
        ii = ii + 1;
        if ii>number_comparisons
            break
        end
        % Simulate observer choice
        if rand()>0.5
            M(pair(1), pair(2)) = M(pair(1), pair(2))+1;
        else
            M(pair(2), pair(1)) = M(pair(2), pair(1))+1;
        end
    end
    
    if trial_i==0
        trial = sum(sum(M))-full_pair_num; % due to initialization 
    else
        trial = sum(M(:));
    end

end

% switch on warning
warning('on','stats:glmfit:IllConditioned')


