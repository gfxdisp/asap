function y = ahess(a,d)

%AHESS  Hessian of scalar field, either audi or function.
%
%   input audi array a:
%   y = AHESS(a)     audi matrix with hessian
%   y = AHESS(a,d) = aeval(ahess(a),0) 
%
%   input function handle f:
%   g = AHESS(f,n) create a function handle g such that g(x1,...,xn) evaluates the 
%   Hessian of f = f(x1,...,xn). Both input and output of g are of class double.
%   g = AHESS(f,var) same, but with an n-vector var of booleans. xi is treated as 
%   variable if var(i) = true and as parameter otherwise.
%
%   Example: <a href="matlab: ahelp(4)">Laplacian of the peaks function and its maximum</a>
%   See also: acurl, adiv, agrad, ajac, alap

if isa(a,'function_handle')
  if isscalar(d)
    d = true(1,d);
  end
  [a1,a2,p] = prepfct(ainit([]),d,2);
  eval(['y = @(' a1 ') ahess(a(' a2 '),0);'])
  return
end

if nargin==2
  y = aeval(ahess(a),d);
  return
end
if numel(a) > 1
  error('Evaluation for scalar functions only.')
else
  y(a.n,a.n) = audi();
  E = eye(a.n);
  for i = 1:a.n
    for j = 1:a.n
      y(i,j) = adiff(a,E(i,:)+E(j,:));
    end
  end
end