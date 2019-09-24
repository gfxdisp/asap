function conditions_order_out=quicksort_rec_sim(conditions_order_in,number_comparisons,q_mat)

    global quickM
    if numel(conditions_order_in) < 2  
        conditions_order_out = conditions_order_in;
        return
    end
    
    pivot = conditions_order_in(end);
    conditions_order_in(end) = [];
 
    ii=0;
    nn=length(conditions_order_in);
    less=[];
    greater=[];
    while (ii<nn) && (sum(quickM(:))<number_comparisons)
        ii=ii+1;
        iia=conditions_order_in(ii);
        
        prob = prob_being_better(q_mat,pivot,iia);
        if rand() < prob
             quickM(pivot,iia)=quickM(pivot,iia)+1;
             less=[less iia];
        else
             quickM(iia,pivot)=quickM(iia,pivot)+1;
             greater=[greater iia];
        end
        
    end
    while (ii<nn)
        ii=ii+1;
        iia=conditions_order_in(ii);
        less=[less iia];
    end
 
    conditions_order_out = [quicksort_rec_sim(less,number_comparisons,q_mat) pivot quicksort_rec_sim(greater,number_comparisons,q_mat)];
 
end
