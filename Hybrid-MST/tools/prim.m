function [mst, cost] = prim(A)

% User supplies adjacency matrix A.  This program uses Prim's algorithm
% to find a minimum spanning tree.  The edges of the minimum spanning
% tree are returned in array mst (of size n-1 by 2), and the total cost
% is returned in variable cost.  The program prints out intermediate 
% results and pauses so that user can see what is happening.  To continue
% after a pause, hit any key.

[n,n] = size(A);                           % The matrix is n by n, where n = # nodes.
%A, n, pause,

if norm(A-A','fro') ~= 0 ,                 % If adjacency matrix is not symmetric,
  disp(' Error:  Adjacency matrix must be symmetric ') % print error message and quit.
  return,
end;

% Start with node 1 and keep track of which nodes are in tree and which are not.

intree = [1];  number_in_tree = 1;  number_of_edges = 0;
notintree = [2:n]';  number_notin_tree = n-1;

in = intree(1:number_in_tree);               % Print which nodes are in tree and which 
out = notintree(1:number_notin_tree); % pause, % are not.

% Iterate until all n nodes are in tree.

while number_in_tree < n,

%   Find the cheapest edge from a node that is in tree to one that is not.

  mincost = Inf;                             % You can actually enter infinity into Matlab.
  for i=1:number_in_tree,               
    for j=1:number_notin_tree,
      ii = intree(i);  jj = notintree(j);
      if A(ii,jj) < mincost, 
        mincost = A(ii,jj); jsave = j; iisave = ii; jjsave = jj;  % Save coords of node.
      end;
    end;
  end;

%    Add this edge and associated node jjsave to tree.  Delete node jsave from list
%    of those not in tree.

  number_of_edges = number_of_edges + 1;      % Increment number of edges in tree.
  mst(number_of_edges,1) = iisave;            % Add this edge to tree.
  mst(number_of_edges,2) = jjsave;
  costs(number_of_edges,1) = mincost;

  number_in_tree = number_in_tree + 1;        % Increment number of nodes that tree connects.
  intree = [intree; jjsave];                  % Add this node to tree.
  for j=jsave+1:number_notin_tree,            % Delete this node from list of those not in tree.
    notintree(j-1) = notintree(j);
  end;
  number_notin_tree = number_notin_tree - 1;  % Decrement number of nodes not in tree.

  in = intree(1:number_in_tree);              % Print which nodes are now in tree and
  out = notintree(1:number_notin_tree); %, pause,% which are not.

end;
%disp(' Edges in minimum spanning tree and their costs: ')
%[mst  costs]                                 % Print out edges in minimum spanning tree.
cost = sum(costs);