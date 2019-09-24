function y = adiv(a,d)

%ADIV Divergence of vector field, either audi or function.
%
%   input audi array a:
%   y = ADIV(a)     audi scalar with divergence of a
%   y = ADIV(a,0) = aeval(adiv(a),0)
%
%   input function handle f:
%   g = ADIV(f,n) create a function handle g such that g(x1,...,xn) evaluates the 
%   divergence of f = f(x1,...,xn). Both input and output of g are of class double.
%   g = ADIV(f,var) same, but with an n-vector var of booleans. xi is treated as 
%   variable if var(i) = true and as parameter otherwise.
%
%   See also: acurl, agrad, ahess, ajac, alap

if isa(a,'function_handle')
  if isscalar(d)
    d = true(1,d);
  end
  [a1,a2,p] = prepfct(ainit([]),d,1);
  eval(['y = @(' a1 ') adiv(a(' a2 '),0);'])
  return
end

if nargin==2
  y = aeval(adiv(a),d);
  return
end
nn = adim(a);
if min(size(a))~=1 || max(size(a))~=nn
  error('Dimension of vector field does not match number of variables.')
else
  E = eye(adim(a));
  y = adiff(a(1),E(1,:));
  for i = 2:adim(a)
    y = y + adiff(a(i),E(i,:));
  end
end