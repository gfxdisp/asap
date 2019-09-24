function n = adim(y)

%ADIM Number of variables.
%   n = ADIM(y) is the number of variables with respect to
%   which derivatives are stored in y.

y1 = y(1);
n = y1.n;