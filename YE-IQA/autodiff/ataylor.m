function p = ataylor(y)

%ATAYLOR Taylor polynomial.
%   p = ATAYLOR(y) is the Taylor polynomial of y = f(x) at x{0}.
%   x is an audi variable with a single point, asize(x) = asize(y) = [1 1].
%   p is a double vector with aord(x)+1 coefficients in descending order.
%   
%   If f has a singularity of type "0/0" at x{0}, p is still computed, but
%   the degree is reduced according to the order of the singularity.
%  
%   Example: <a href="matlab: ahelp(5)">Sinc function and rule of l'Hospital</a>
%   See also: polyval, aode

p = intern_taylor(y);