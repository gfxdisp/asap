function s = asize(y)

%ASIZE Size of data stored in an audi variable.
%   s = ASIZE(y) is the size of the arrays containing function
%   values and derivatives of y.

s = y(1).s;