function [prob] = prob_being_better(q,id1,id2)

    if size(q,1)>1 && size(q,2)>1
        prob = (q(id1,id2)+1)/(q(id1,id2)+q(id2,id1)+2);
    else
        prob = normcdf( q(id1)-q(id2), 0, 1.4865);
    end
end