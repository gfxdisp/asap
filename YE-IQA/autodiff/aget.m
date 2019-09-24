function b = aget(a,varargin)

%AGET Select parts of audi data.
%   b = AGET(a,i) or b = AGET(a,i,j) read the function values and
%   derivatives with indices i or i,j and store it in the audi variable b.
%
%   See also: aset

b = a;
if numel(a) > 1
  for i = 1:numel(a)
    b(i) = aget(a(i),varargin{:});
  end
else
  if length(varargin)==1
    i = varargin{1};
    if islogical(i)
      i = find(i);
    end
    b.s = size(i);
    b.c = a.c(i(:),:);
  else
    i = varargin{1};
    j = varargin{2};
    if islogical(i)
      i = find(i);
    end
    if islogical(j)
      j = find(j);
    end
    I = ones(a.s);
    I(i,j) = 0;
    b.s = [numel(i) numel(j)];
    b.c = a.c(~I,:);
  end
end