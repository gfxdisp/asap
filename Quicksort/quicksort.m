function conditions_order_out=quicksort(conditions_order_in,number_of_comparisons)

    global quickM

    if numel(conditions_order_in) < 2  
        conditions_order_out = conditions_order_in;
        return
    end
    
    pivot = conditions_order_in(end);
    conditions_order_in(end) = [];
 
    ii=0;
    number_conditions=length(conditions_order_in);
    less=[];
    greater=[];
    
    while (ii<number_conditions) && (sum(quickM(:))<number_of_comparisons)
        ii=ii+1;
        iia=conditions_order_in(ii);
        
        % 
        % Simulate observer choice, use the vote as a pivot in sorting
        % Change from here 
        if rand() < 0.5
             quickM(pivot,iia)=quickM(pivot,iia)+1;
             less=[less iia];
        else
             quickM(iia,pivot)=quickM(iia,pivot)+1;
             greater=[greater iia];
        end
        % Till here 
        
    end
    while (ii<number_conditions)
        ii=ii+1;
        iia=conditions_order_in(ii);
        less=[less iia];
    end
    conditions_order_out = [quicksort(less,number_of_comparisons) pivot quicksort(greater,number_of_comparisons)];
 
end
