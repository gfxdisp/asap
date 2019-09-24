%% Automatic Differentiation with the *AutoDiff* toolbox
% 
% by 
% <http://www3.mathematik.tu-darmstadt.de/hp/geometrie-und-approximation/reif-ulrich/home.html 
% Ulrich Reif>
%
% Version 3.1 - March 31, 2017
%
% The *AutoDiff toolbox* realizes automatic differentiation of functions
% by means of operator overloading. Without relying on the symbolic toolbox,
% it provides 
%
% * Ordinary and partial derivatives of arbitrary order
% * Differential operators, like gradient, Laplacian, or curl
% * Creation of stand-alone functions, as requested for optimization or ODE tools
% * Taylor expansions of explicit expressions and solutions of ODEs
% * Curvature computation for curves and surfaces
% * Resolution of singularities of type "0/0" by l'Hospital's rule
%
% Automatic differentiation is useful whenever derivatives of
% functional expressions are needed and calculus by hand is too tedious. 
% Finite differences and symbolic computations are alternatives with
% pros and cons known to anybody working with them.
%
% The _audi_ class, which is the central part of the toolbox, admits to 
% evaluate derivatives of functions which
% are defined in the usual way by composing arithmetic operators and
% elementary functions. During evaluation, derivatives are computed in 
% the background without further user interaction if the input variables 
% are of class _audi_. In principle, automatic differentiation yields 
% exact results, but of course, accuracy is subject to rounding errors, 
% as for any other numerical computation.
%
% This release is based on a complete redesign of the _audi_ class.
% While being fully downward compatible, it offers
% 
% * a significant speedup
% * new possibilities for dealing with matrices, like |eig,chol,pinv,...|
% * many more overloaded methods
% * singleton expansion for basic arithmetic operations.
%
% This version of *AutoDiff* requests Matlab Release R2016b or higher.
% Version 2.2, which runs with earlier releases of Matlab, is still 
% available, but will not be developed any further.

%% Initialization and evaluation
% When initializing arguments of functions as _audi_ variables, subsequent
% results of computations do not only contain functions values, as usual,
% but also derivatives up to some given order |k|. 
%
% *Univariate case:* 
% Let |f| be a function depending on one variable.
% To evaluate |f| and all its derivatives up to order |k|
% at some points |v|, stored as a _double_ array, type
%
%       x = ainit(v,k)
%       y = f(x)
%
% Curly braces are used to get the values |yd| of the |d|-th derivative
% of |f| at the points |v|:
%
%       yd = y{d}
%
% The size of the double array |yd| is
%
%       asize(y) = size(yd) = size(v)

%%

% Example
x = ainit(1:5,3);      % initialize audi variable, up to 3rd derivative covered
y = x.^3;              % compute function f(x) = x^3
y{2}                   % values of 2nd derivative, f''(x) = 6x

%%
% *Multivariate case:* If |f| depends on |n| variables and shall be evaluated
% at the some points |v1,...,vn|, stored as _double_ arrays of equal size, type
%
%      [x1,...,xn] = ainit(v1,...,vn,k)
%      y = f(x1,...,xn)
%
% Again, curly braces are used to get the values |yd| of the partial derivative of order
% |d1,...,dn| at the points |v1,...,vn|:
%
%       yd = y{d1,...,dn}
%
% The sum |d1+...+dn| must not exceed
%
%       aord(y) = k
%
% Curly braces can be replaced by the command |aeval|,
%
%       y{d1,...,dn} = aeval(y,d1,...,dn)
%
% Use the command |amatinit| to initialize matrices containing a larger number of 
% independent _audi_ variables.
%% 

% Example
[u,v] = ainit([1 -1 2],[6 2 4],4);  % initialize two independent audi variables
w = v/u;                            % compute function w(u,v) = v/u
w{2,1}                              % values of partial derivative w_uuv

z = amatinit([3 7 2 4],1);          % vector with 4 independent audi variables
aeval(norm(eye(4)+z*z'),[0,0,1,0])  % derivative of 2-norm of Id+z*z' wrt z(3)

%% 
% 
% _Example:_ 
% <example1.html Functions, derivatives, and Newton's method> 
%
% _See also:_ 
% <matlab:doc('ainit') ainit>,
% <matlab:doc('amatinit') amatinit>,
% <matlab:doc('aeval') aeval>,
% <matlab:doc('asize') asize>,
% <matlab:doc('aord') aord>

%% Arrays
% _audi_ variables can be concatenated to form arrays in the familiar way.
% When combining _audi_ variables and _doubles_, the latter are treated as
% constants.
% Also operations like |repmat| or |transpose| are supported. _audi_ arrays
% can be multiplied by other arrays of class _double_ or _audi_ in the sense
% of matrix and elementwise multiplication using the operators |*| and |.*|,
% respectively.
%
% Parentheses select subarrays of _audi_ arrays. By contrast, to define the 
% _audi_ array |b| which contains only parts of the data of the given 
% _audi_ array |a|, use
%
%       b = aget(a,ref) 

%
% where |ref| is a single index or a list of indices addressing parts 
% of the data arrays. Note that |size(b)=size(a)|, while |asize(b)| and
% |asize(a)| may be different.
% Parts of _audi_ data are replaced by
%
%       a = aset(a,b,ref) 
%
% Evaluation of arrays by curly braces follows special formatting rules, as
% documented for |aeval|.
%
% _Example:_ 
% <example2.html Curvature and torsion of a space curve> 
%
% _See also:_ 
% <matlab:doc('aeval') aeval>,
% <matlab:doc('aget') aget>,
% <matlab:doc('aset') aset>

%% Taylor expansion of explicit functions and ODEs
% *Explicit case:* Let |f| be a scalar function of one variable. To determine
% the coefficients |p| of the Taylor polynomial of degree |k| at the point |x0|, 
% type
%
%       x = ainit(x0,k)
%       p = ataylor(f(x))

%%

% Example
x = ainit(1,4);
p = ataylor(log(x))

%%
% *Solution of ODE:* Consider the initial value problem
%
%       y'(x) = f(x,y)
%       y(x0) = y0
%
% To obtain the coefficients |p| of the degree |k| Taylor polynomial of the sought
% function |y| at |x0|, type
%
%       p = aode(f,x0,y0,k)
%
% where |f = @(x,y)...| is an anonymous function. For systems of ODEs,
% the Taylor polynomials |p1,...,pm| of the sought functions |y1,...,ym|
% can be computed by
%
%       [p1,...,pm] = aode(f,x0,y0,k)
%
% where |f = @(x,y1,...,ym)[...;...;...]| is an anonymous function returning
% a column vector with |m| entries.
%
% It is possible to pass additional parameters to the routine,
%
%       [p1,...,pm] = aode(f,x0,y0,k,c1,c2,...)
%
% _Example:_ 
% <example3.html The mathematical pendulum> 
%
% _See also:_ 
% <matlab:doc('ataylor') ataylor>,
% <matlab:doc('aode') aode>

%% Differential operators
%
% Differential operators map _audi_ variables to new _audi_ variables of lower order.
% 
% The general differentiation operator is
%
%       yd = adiff(y,d1,...,dn)
%
% The order of the _audi_ variable |yd| is
%
%       aord(yd) = aord(y) - (d1+...+dn)
%
% The _audi_ class also provides various other predefined differential operators:
%
% * scalar-valued: |alap| (Laplacian), |adiv| (divergence)
% * vector-valued: |agrad| (gradient), |acurl| (curl)
% * matrix-valued: |ajac| (Jacobian), |ahess| (Hessian)
% 
% There exists a shortcut to evaluate the outcome of any of these operators:
%
%      a*(y,0) = aeval(a*(y),0,...,0)
% 
% 
% _Example:_ 
% <example4.html  Laplacian of the peaks function and its maximum> 
%
% _See also:_ 
% <matlab:doc('acurl') acurl>,
% <matlab:doc('adiff') adiff>,
% <matlab:doc('adiv') adiv>,
% <matlab:doc('agrad') agrad>,
% <matlab:doc('ahess') ahess>,
% <matlab:doc('ajac') ajac>,
% <matlab:doc('alap') alap>

%% Stand-alone functions
%
% Sometimes, it is convenient or even necessary to have available a stand-alone 
% function for evaluating gradients or other differential properties of a given 
% function. In particular, if you want to 
%
% * pass derivatives to routines like |integral| or |fzero|
% * enhance the performance of optimization routines like |fmincon| by supplying
%   gradients or Hessians
% * furnish solvers for stiff ODEs like |ode15s| by supplying Jacobians 
%
% such functionality is requested. To this end, all differential operators described 
% in the preceding section are able to create function handles from function handles. 
% Using this feature, initialization of _audi_ variables becomes obsolete at the
% cost of slightly increased computation times.

%%

% Example
f_uvv = adiff(@peaks,[1 2]);    % 2nd argument = order of differentiation
f_uvv(-2:2,zeros(1,5))          % evaluate partial derivative

L = alap(@peaks,2);             % 2nd argument = number of variables
L(-2:2,zeros(1,5))              % evaluate Laplacian

%%
% It is important to note that the arguments to the functions |f_uvv| and |L|
% are standard double arrays, and not _audi_ variables.

%%
% When solving a stiff ODE |y' = f(t,y),| the Jacobian |J| of |f| with respect to
% |y| may be supplied to solvers like |ode15s| to enhance performance. The
% according function handle can be created as follows:
%
%       g = @(t,y1,...,yn) f(t,[y1;...;yn]);   % function with scalar arguments
%       j = ajac(g,[0 1...1]);                 % Jacobian, 1st argument is parameter
%       J = @(t,y) j(t,y(1),...,y(n));         % Jacobian with vector argument
%
% According constructions apply for constrained and unconstrained optimization.
%
% _See also:_ 
% <matlab:doc('acurl') acurl>,
% <matlab:doc('adiff') adiff>,
% <matlab:doc('adiv') adiv>,
% <matlab:doc('agrad') agrad>,
% <matlab:doc('ahess') ahess>,
% <matlab:doc('ajac') ajac>,
% <matlab:doc('alap') alap>,
% <matlab:doc('ode15s') ode15s>,
% <matlab:doc('fmincon') fmincon>


%% Rule of l'Hospital
%
% Singularities of type |"0/0"| are resolved automatically using the 
% rule of l'Hospital as long as the order of the singularity does not 
% exceed the order of the _audi_ variable. Information on higher order 
% derivatives can get lost and is then replaced by |NaN|.
%
% _Example:_ 
% <example5.html Sinc function and rule of l'Hospital> 
%
% See also: 
% <matlab:doc('ataylor') ataylor>

%% Matrix computations
%
% The _audi_ class provides an extended set of overloaded methods for dealing
% with matrices, including the solution of linear systems and factorization. 
% However, not all options provided by the corresponding built-in routines are 
% supported. Further, it is important to note that all routines assume generic 
% input. For instance, the null space of a square matrix is always empty, even 
% if all elements are zero.

%%
% Example:

x = amatinit([-2 1 2 -1],2);      % vector with four independent audi variables
A = x'*x - 2*diag(x);             % 4x4 matrix
ahess(max(abs(eig(inv(A)))),0)    % Hessian of the spectral radius of the inverse 

%% 
%
% _Example:_ 
% <example6.html Gaussian and mean curvature of a Klein bottle> 

%% Curves and Surfaces
%
% Differential geometric properties of parametrized curves and surfaces can 
% be computed conveniently using the programs |acurve| and |asurf|, respectively.
%
% Curvature as well as tangent and normal vector of a planar curve |C| are given by
%
%      [kap,T,N] = acurve(C)
%
% If |C| is a space curve, also torsion and the binormal vector are provided,
%
%      [kap,tau,T,N,B] = acurve(C)
%
% For a surface |S|, the Gaussian and mean curvature, the normal vector, and the
% first two fundamental forms are computed by
%
%      [K,H,N,G,B] = asurf(S)
%
% _Example:_ 
% <example7.html Curvature revisited> 
%
% See also: 
% <matlab:doc('acurve') acurve>,
% <matlab:doc('asurf') asurf>

%% Overloaded operators and functions
%
% * unary and binary arithmetic operators: all
% * relational and logical operators: all but |&&| and | || |
% * elementary functions: |exp,sin,cos,tan,sinh,cosh,tanh,sqrt| and their inverses
% * vectors and matrices: |det,trace,inv,pinv,norm,eig,svd,lu,qr,chol,null,orth|
%   |sum,cumsum,prod,cumprod,diff,diag,tril,triu,cross,dot,expm,logm,sqrtm|
% * polynomials: |polyval,polyvalm,polyfit,polyder,polyint,conv,deconv|
% * classification: |isreal,isdiag,istriu,istril,issymmetric,ishermitean,isfinite|
% * complex: |abs,real,imag,conj,angle|
% * miscellaneous: |min,max,sign,fft,ifft,all,any,bsxfun,atan2|
%
% To overload any other function |fct| defined for doubles, insert the lines
%
%       function y = fct(x)
%         df = @(x) ...;
%         y  = adval(df,x);
%       end
%
% in the |methods| section of the class definition |audi.m|. The anonymous 
% function |df| must return the derivative of |fct| at |x|. 
%
% If |fct| depends on additional parameters, say |fct(a,b,x,c)|, then use the syntax
%
%       function y = fct(a,b,x,c)
%         df = @(a,b,x,c) ...;
%         y  = adval(df,a,b,x,c,3);
%       end
%
% to indicate that |df| is the partial derivative of |fct| with respect
% to the third argument.
%
% _Example:_ 
% See the methods |sin| and |polyval| in the class definition |audi.m|.


%% Limitations
%
% * The set of overloaded methods is not complete, but should
%   cover a substantial class of applications. Additional functions are readily 
%   overloaded when needed.
% * The implementation aims at generality and ease of use at the cost of a
%   certain slow-down, in particular when requesting high order derivatives
%   for functions with many unknowns.
% * Currently, there is no mechanism to overload built-in functions depending
%   on several variables, like |ellipj|. In such cases, one argument can be
%   selected as variable, while all others are treated as parameters.
% * Function input is rarely checked for validity. So error messages
%   may be confusing when evoked by subordinate functions.
% * Matrix routines like |qr| or |svd| do not cover all variants of their
%   built-in counterparts. Also complex matrices are not fully supported.

%% Release notes
%
% * AutoDiff 1.0, 05/01/16: 
%      launch
% * AutoDiff 2.0, 05/09/16: 
%      curly braces for arrays;
%      stand-alone functions;
%      |aode| fixed; 
%      min/max yield derivatives from the right;
%      new overloaded functions |atan2|, |real|, |imag|, |conj|, |angle|, |isfinite|, 
%      |isreal|
% * AutoDiff 2.1, 10/30/16:
%      initialization of vectors of audi variables with |avecinit|;
%      negative integer exponents for |mpower|;
%      concatenation of _audi_ and _double_ arrrays;
%      new overloaded functions |bsxfun|, |fft|, |ifft|;
%      arbitrary variable names in |aode|;
%      bug fix in |acurve|;
% * AutoDiff 3.0, 03/02/17:
%      redesign of the class; 
%      significant speedup;
%      many more overloaded methods, in particular for matrix computations;
%      singleton expansion for arithmetic operators;
%      minor bug fixes.
% * AutoDiff 3.1, 03/31/17:
%      speedup for stand-alone functions.