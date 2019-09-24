function varargout = acurve(C)

%ACURVE Differential geometric properties of a curve.
%   [kap,T,N] = ACURVE(C) determines some differential geometric
%   properties of the planar curve C, given as an audi column
%   vector containing the two coordinate functions. The order 
%   of S must be at least 2. The output arguments are audis with
%
%      kap - signed curvature
%      T   - tangent vector
%      N   - normal vector
%
%   [kap,tau,T,N,B] = ACURVE(C) determines some differential geometric
%   properties of the space curve C, given as an audi column
%   vector containing the three coordinate functions. The order
%   of S must be at least 3. The output arguments are audis with
%
%      kap - curvature
%      tau - torsion
%      T   - tangent vector
%      N   - normal vector
%      B   - binormal vector
%
%   Example: <a href="matlab: ahelp(7)">Curvature revisited</a>
%   See also: asurf

if numel(C) == 2
  c1  = adiff(C,1);
  c2  = adiff(C,2);
  T   = c1/norm(c1);
  N   = [-T(2);T(1)];
  kap = (c1(1)*c2(2) - c1(2)*c2(1))/norm(c1).^3;
  varargout = {kap,T,N};
elseif numel(C) == 3
  c1  = adiff(C,1);
  c2  = adiff(C,2);
  b   = cross(c1,c2);
  kap = sqrt(sum(b.^2)./sum(c1.^2).^3);   
  tau = dot(b,adiff(C,3))./sum(b.^2);
  T   = c1/norm(c1);
  B   = cross(T,c2);
  B   = B/norm(B);
  N   = cross(B,T);
  varargout = {kap,tau,T,N,B};
end