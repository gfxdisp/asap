function y = adiff(a,varargin)

%ADIFF Partial derivative of audi array or function.
%
%   input audi array a:
%   y = ADIFF(a,0)            do nothing
%   y = ADIFF(a, d1,...,dn)   partial derivative of order d1,...,dn
%   y = ADIFF(a,[d1,...,dn])  same with order vector
%
%   input function handle f:
%   g = ADIFF(f,[d1,...,dn])  create a function handle g such that g(x1,...,xn) 
%   evaluates the partial derivative of order d1,...,dn. Both input and output 
%   of g are of class double.
%
%   Example: <a href="matlab: ahelp(6)">Gaussian and mean curvature of a Klein bottle</a>
%   See also: acurl, adiv, agrad, ahess, ajac, alap

if isa(a,'function_handle')
  d = cell2mat(varargin);
  [a1,a2,p] = prepfct(ainit([]),d);
  eval(['y = @(' a1 ') aeval(adiff(a(' a2 '),d));'])
  return
end

if numel(a)>1
  y = a;
  for i = 1:numel(a)
    y(i) = adiff(a(i),varargin{:});
  end
  return
end
if nargin == 1
  d = zeros(1,aord(a));
else
  d = cell2mat(varargin);
end
y = a;
y.k = a.k - sum(d);
y.c = a.c(:,1:nchoosek(y.n+y.k,y.n));
audi([],0,adim(a),aord(a)-sum(d));
for i = 1:size(y.c,2)
  y.c(:,i) = a.c(:,idx(a,sub(y,i)+d));
end