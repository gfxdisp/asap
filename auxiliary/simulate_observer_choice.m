function [MG] = simulate_observer_choice(q,id1,id2,MG)

    prob = prob_being_better(q,id1,id2);
    if rand() < prob
        if size(MG,1)==size(MG,2) && size(MG,1)==numel(q)
            MG(id1,id2)=MG(id1,id2)+1;
        else
            MG = [MG; id1,id2];
        end
    else
        if size(MG,1)==size(MG,2) && size(MG,1)==numel(q)
            MG(id2,id1)=MG(id2,id1)+1;
        else
            MG = [MG; id2,id1];
        end
    end
end
