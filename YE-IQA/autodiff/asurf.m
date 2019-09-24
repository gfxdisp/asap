function [K,H,N,G,B] = asurf(S)

%ASURF Differential geometric properties of a surface.
%   [K,H,N,G,B] = ASURF(S) determines some differential geometric
%   properties of the parametrized surface S, given as an audi column
%   vector containing the three coordinate functions. The order of S 
%   must be at least 2. The output arguments are audis with
%
%      K - Gaussian curvature
%      H - mean curvature
%      N - normal vector
%      G - 1st fundamental form
%      B - 2nd fundamental form
%
%   Example: <a href="matlab: ahelp(7)">Curvature revisited</a>
%   See also: acurve

J = ajac(S);                              
N = cross(J(:,1),J(:,2));                
N = N/norm(N);                           
G = J'*J;                               
B =-J'*ajac(N);                                 
K = det(B)/det(G);                              
H = trace(G\B)/2;