

number_of_conditions = 10;
number_of_comparisons = 400;

C_init = zeros (number_of_conditions, number_of_conditions);
ii = 0;
while ii<number_of_comparisons

    [pairs_to_compare] = ASAP(C_init, 'mst');
    
    for jj = 1:size(pairs_to_compare,1)
        
        
        pair = pairs_to_compare(jj,:);
        disp(['Round: ',num2str(ii),', observer compares: ',num2str(pair(1)),' ',num2str(pair(2))])
        
        ii = ii + 1;
        if ii > number_of_comparisons
            break
        end
        
        % Simulate observer choice
        if rand()>0.5
            C_init(pair(1), pair(2)) = C_init(pair(1), pair(2))+1;
        else
            C_init(pair(2), pair(1)) = C_init(pair(2), pair(1))+1;
        end
    end
end