function y = alap(a,d)

%ALAP Laplacian of scalar field, either audi or function.
%
%   input audi a:
%   y = ALAP(a)     audi with Laplacian
%   y = ALAP(a,0) = aeval(alap(a),0)
%
%   input function handle f:
%   g = ALAP(f,n) create a function handle g such that g(x1,...,xn) evaluates the 
%   Laplacian of f = f(x1,...,xn). Both input and output of g are of class double.
%   g = ALAP(f,var) same, but with an n-vector var of booleans. xi is treated as 
%   variable if var(i) = true and as parameter otherwise.
%
%   Example: <a href="matlab: ahelp(4)">Laplacian of the peaks function and its maximum</a>
%   See also: acurl, adiv, agrad, ahess, ajac

if isa(a,'function_handle')
  if isscalar(d)
    d = true(1,d);
  end
  [a1,a2,p] = prepfct(ainit([]),d,2);
  eval(['y = @(' a1 ') alap(a(' a2 '),0);'])
  return
end

if nargin==2
  y = aeval(alap(a),d);
  return
end
if numel(a) > 1
  y = a;
  for i = 1:numel(a)
    y(i) = alap(a(i));
  end
  return
else
  E = eye(a.n);
  y = adiff(a,2*E(1,:));
  for i = 2:a.n
    y = y + adiff(a,2*E(i,:));
  end
end