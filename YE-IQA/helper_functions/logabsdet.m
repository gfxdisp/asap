function v = logabsdet(A, op)
%LOGABSDET Computation of logarithm of determinant of a matrix
%
%   v = logabsdet(A);
%       computes the logarithm of the absolute value of the determinant of A. 
%
%       Here, A should be a square matrix of double or single class.
%       If A is singular, it will returns -inf.
%
%       Theoretically, this function should be functionally 
%       equivalent to log(abs(det(A))). However, it avoids the 
%       overflow/underflow problems that are likely to 
%       happen when applying det to large matrices.
%
%       The key idea is based on the mathematical fact that
%       the determinant of a triangular matrix equals the
%       product of its diagonal elements. Hence, the matrix's
%       log-abs-determinant is equal to the sum of the absolute values of
%       their logarithms. By keeping all computations in log-scale, the
%       problem of underflow/overflow caused by product of 
%       many numbers can be effectively circumvented.
%
%       The implementation is based on LU factorization.
%
%   v = logdet(A, 'chol');
%       If A is positive definite, you can tell the function 
%       to use Cholesky factorization to accomplish the task 
%       using this syntax, which is substantially more efficient
%       for positive definite matrices. 
%
%   Remarks
%   -------
%       logarithm of the absolute value of the determinant of a matrix
%       widely occurs in the context of multivariate statistics. This
%       function might be useful there, especially in a high-dimensional
%       space.       
%
%       Theoretially, LU, QR can both do the job. However, LU 
%       factorization is substantially faster. So, for generic
%       matrices, LU factorization is adopted. 
%
%       For positive definite matrices, such as covariance matrices,
%       Cholesky factorization is typically more efficient. And it
%       is STRONGLY RECOMMENDED that you use the chol (2nd syntax above) 
%       when you are sure that you are dealing with a positive definite
%       matrix.
%
%   Examples
%   --------
%       % compute the log-abs-determinant of a generic matrix
%       A = rand(1000);
%       v = logabsdet(A);
%
%       % compute the log-abs-determinant of a positive-definite matrix
%       A = rand(1000);
%       C = A * A';     % this makes C positive definite
%       v = logabsdet(C, 'chol');
%

%   Copyright 2008, Dahua Lin, MIT
%   Email: dhlin@mit.edu
%   Modified by Clifford Wilmot, University of Cambridge, 2018
%
%   This file can be freely modified or distributed for any kind of 
%   purposes.
%

%% argument checking

assert(isfloat(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
    'logdet:invalidarg', ...
    'A should be a square matrix of double or single class.');

if nargin < 2
    use_chol = 0;
else
    assert(strcmpi(op, 'chol'), ...
        'logdet:invalidarg', ...
        'The second argument can only be a string ''chol'' if it is specified.');
    use_chol = 1;
end

%% computation

if use_chol
    v = 2 * sum(log(abs(diag(chol(A)))));
else
    [L, U, P] = lu(A);
    du = diag(U);
    c = det(P) * prod(sign(du));
    v = log(abs(c)) + sum(log(abs(du)));
end