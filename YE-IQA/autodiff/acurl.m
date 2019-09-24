function y = acurl(a,d)

%ACURL Curl of a 3d vector field, either audi or function.
%
%   input audi array a:
%   y = ACURL(a)     audi 3d column vector with curl of a
%   y = ACURL(a,0) = aeval(curl(a),0)
%
%   input function handle f:
%   g = ACURL(f)     create a function handle g such that g(x,y,z) 
%   evaluates the curl of f = f(x,y,z). Both input and output of g are of 
%   class double.
%   g = ACURL(f,var) same, but for f = f(x1,...,xn) and an n-vector var of 
%   booleans. xi is treated as variable if var(i) = true and as parameter
%   otherwise.
%
%   See also: adiv, agrad, ahess, ajac, alap

if isa(a,'function_handle')
  if nargin == 1
    d = true(1,3);
  end
  [a1,a2,p] = prepfct(ainit([]),d,1);
  eval(['y = @(' a1 ') acurl(a(' a2 '),0);'])
  return
end


if nargin==2
  y = aeval(acurl(a),d);
  return
end
if ~isequal(size(a),[3,1])
  error('Size of vector field must be [3 1].')
else
  J = ajac(a);
  y = [J(3,2)-J(2,3); J(1,3)-J(3,1); J(2,1)-J(1,2)];
end