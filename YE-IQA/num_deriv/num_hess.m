% num_hess.m
%
% H = num_hess(func, X, h)
%
% Function to compute the numerical Hessian of an arbitrary objective
% function.
%
% Inputs:
%	-> func: Function handle for which numerical Hessian is to be
%		obtained.
%	-> X: Point of interest about which Hessian is to be obtained.
%	-> h: Tolerance for differentiation.
%
% Output:
%	-> H: Numerical Hessian of function func (square symmetric
%		matrix of size n=length(X)).
%
% Created by: Brendan C. Wood
% Created on: February 14, 2011
%
% Copyright (c) 2011, Brendan C. Wood <b.wood@unb.ca>


function H = num_hess(func, X, h)

H = zeros(length(X));

% for each dimension of objective function
for i=1:length(X)
    % derivative at first point (left)
    x1 = X;
    x1(i) = X(i) - h;
    df1 = num_grad(func, x1, h);
    
    % derivative at second point (right)
    x2 = X;
    x2(i) = X(i) + h;
    df2 = num_grad(func, x2, h);
    
    % differentiate between the two derivatives
    d2f = (df2-df1) / (2*h);
    
    % assign as row i of Hessian
    H(i,:) = d2f';
end