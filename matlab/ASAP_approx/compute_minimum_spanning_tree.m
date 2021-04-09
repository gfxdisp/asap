function [pairs_to_compare] = compute_minimum_spanning_tree(inf_mat)
    % Asap allows for batch mode, for which we return N-1 pairs of conditions for comparisons
    % forming the minimum spanning tree of the 1/information_gain matrix
    inf_mat = inf_mat+inf_mat';
    inf_mat = 1./inf_mat;
    inf_mat(inf_mat<0) = Inf;
    GrMST = graph(inf_mat);

    t = minspantree(GrMST);

    t_edges = sortrows(t.Edges,2);
    pairs_to_compare = t_edges.EndNodes;
end
