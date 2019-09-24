function varargout = ainit(varargin)

%AINIT Initialize audi variable.
%   x = AINIT(v) is an audi variable with value v and derivative 1 
%   to be used as argument of univariate functions, as in y = f(x). 
%   v is a double array of arbitrary dimension.
%
%   x = AINIT(v,k) is an audi variable with derivatives up to order k, 
%   where all derivatives of order 2 or higher are set to 0. 
%
%   [x1,...,xn] = AINIT(v1,...,vn,[k=1]) is a set of audi variables to
%   be used as arguments of functions depending on n variables, as in 
%   y = f(x1,...,xn).
%
%   x = AINIT(v,n,d,k) is an audi variable of dimension n and order k which
%   is all zero, except for x{d} = v.
%
%   Example: <a href="matlab: ahelp(1)">Functions, derivatives, and Newton's methods</a>
%   See also: avecinit, aeval, adim, aord, asize

if nargout==1 && nargin==4
  x = 0*audi(0*varargin{1},1,varargin{2},varargin{4});
  x.c{idx(x,varargin{3})} = varargin{1};
  varargout{1} = x;
  return
end

n = nargout;
if nargin == n
  k = 1;
else
  k = varargin{end};
end
varargout = cell(1,n);
for i = 1:n
  if isa(varargin{i},'audi')
    varargin{i} = aeval(varargin{i});
  end
  if ~isequal(size(varargin{1}),size(varargin{i}))
    error('All variables must have equal size.')
  end
  varargout{i} = audi(varargin{i},i,n,k);
end