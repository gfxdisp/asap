function x = avecinit(v,k)

%AVECINIT Initialize audi vector.
%   AVECINIT is obsolete and will be removed in a future 
%   release. Use amatinit instead.
%
%   x = AVECINIT(v) for a double array v is an audi array 
%   of n=numel(v) independent variables of order k with 
%   size(x) = size(v) and values x{0} = v.
%
%   See also: amatinit, ainit

warning('avecinit is obsolete and will be removed in a future release. Use amatinit instead')
if nargin==1
  k = 1;
end
for i = numel(v):-1:1
  x(i,1) = audi(v(i),i,numel(v),k);
end
x = reshape(x,size(v));