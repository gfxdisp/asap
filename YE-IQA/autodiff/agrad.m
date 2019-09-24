function y = agrad(a,d)

%AGRAD Gradient of scalar field, either audi or function.
%
%   input audi array a:
%   y = AGRAD(a)     audi column vector with gradient
%   y = AGRAD(a,0) = aeval(agrad(a),0)
%
%   input function handle f:
%   g = AGRAD(f,n) create a function handle g such that g(x1,...,xn) evaluates the 
%   gradient of f = f(x1,...,xn). Both input and output of g are of class double.
%   g = AGRAD(f,var) same, but with an n-vector var of booleans. xi is treated as 
%   variable if var(i) = true and as parameter otherwise.
%
%   Example: <a href="matlab: ahelp(4)">Laplacian of the peaks function and its maximum</a>
%   See also: acurl, adiv, ahess, ajac, alap

if isa(a,'function_handle')
  if isscalar(d)
    d = true(1,d);
  end
  [a1,a2,p] = prepfct(ainit([]),d,1);
  eval(['y = @(' a1 ') agrad(a(' a2 '),0);'])
  return
end

if nargin==2
  y = aeval(agrad(a),d);
  return
end
if numel(a) > 1
  error('Evaluation for scalar functions only.')
else
  E = eye(a.n);
  y(a.n,1) = audi();
  for i = 1:a.n
    y(i) = adiff(a,E(i,:));
  end
end