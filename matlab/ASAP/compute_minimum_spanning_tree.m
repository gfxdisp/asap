function [pairs_to_compare] = compute_minimum_spanning_tree(inf_mat)
    inf_mat = inf_mat+inf_mat';
    inf_mat = 1./inf_mat;
    inf_mat(inf_mat<0) = Inf;
    GrMST = graph(inf_mat);

    t = minspantree(GrMST);

    t_edges = sortrows(t.Edges,2);
    pairs_to_compare = t_edges.EndNodes;
end