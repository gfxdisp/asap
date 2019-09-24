function [M] = quicksort_sim(q,N,number_comparisons)
    
    global quickM
    quickM = zeros(N,N);
    
    conditions_order=randperm(N);
    
    while sum(quickM(:))<quickcap
        conditions_order=quicksort_rec_sim(conditions_order,number_comparisons,q);
    end
    M = quickM;
end

