function a = aset(a,b,varargin)

%ASET Overwrite data.
%   a = ASET(a,b,i) or a = ASET(a,b,i,j) replace the function values and
%   derivatives with indices i or i,j by those in the audi variable b.
%
%   See also: aget

if numel(varargin)==1
  a.c(varargin{1},:) = b.c;
else
  i = varargin{1};
  j = varargin{2};
  I = ones(a.s);
  I(i,j) = 0;
  a.c(~I,:) = b.c;
end