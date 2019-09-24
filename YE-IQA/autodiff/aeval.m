function v = aeval(y,varargin)

%AEVAL Evaluate audi array.
%   v = AEVAL(y) or v = AEVAL(y,0) yield a double array with the function values 
%   stored in the audi array y. 
%   v = AEVAL(y,d1,...,dn) or v = AEVAL(y,[d1,...,dn]) yield the values of the
%   derivative of order d1,...,dn stored in y.
%   AEVAL can be replaced by curly braces, y{d1,...,dn} = AEVAL(y,d1,...,dn).
%
%   The following formatting rules apply:
%   * If the values are scalars, asize = [1 1], then each element of the array 
%     y corresponds to one element of v, size(v) = size(y).
%   * If y is a single audi variable, size(y) = [1 1], then v is the array of
%     values stored in y with size(v) = asize(y).
%   * If y is a row vector of m audi variables with always N = prod(asize(y)) 
%     values, then y is an Nxm matrix. Linear indexing applies.
%   * If y is a column vector of m audi variable with always N = prod(asize(y))
%     values, then y is a mxN matrix. Linear indexing applies.
%   * Otherwise, if y is an array of dimension D, then v is an array of dimension 
%     D+1, where the last coordinate contains N function values obtained by linear 
%     indexing.
%
%   See also: ainit, asize

if ~isa(y,'audi')
  v = y;
  return
end
b = y(1);
if nargin==1 || isequal(varargin{:},0)
  d = zeros(1,b.n);
else
  if isempty(varargin{end})
    d = zeros(1,y(1).n);
    for i = 1:numel(varargin)-1
      for j = 1:numel(varargin{i})
        d(varargin{i}(j)) = d(varargin{i}(j)) + 1;
      end
    end
  else
    d = cell2mat(varargin);
  end
  if sum(d)>aord(y) || ~isequal(size(d),[1 adim(y)])
    error('Requested derivative not available.')
  end
end
I = idx(b,d);
if numel(y) == 1
  v = reshape(y.c(:,I),y.s);
elseif size(y,1) == 1
  v = (aeval(y',d))';
else
  s = size(y);
  if s(end)==1
    s = s(1:end-1);
  end
  v = zeros([numel(b.c(:,1)) s]);
  for i = 1:numel(y)
    b = y(i);
    z = b.c(:,I);
    v(:,i) = z(:);
  end
  v = permute(v,[2:numel(s)+1 1]);
end