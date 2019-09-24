% This function uses automatic differentiation to compute the Hessian of an
% anonymous function f(x) at the point x.
function h=quick_ahess(f,x)
a=num2cell(x);
% Must assign ainit variables at same time.
adims=size(a);
N=adims(2);
a2=cell(1,N);
[a2{:}]=ainit(a{:},2);

% Define auto-diff version of f(x) for the values given in a.
g=f([a2{:}]);

% Evaluates the Hessian at the values given in a (the 0 means that we
% take 0th order derivatives of the Hessian (i.e. the Hessian itself)).
% This is equivalent to:
% h=ahess(g)
% h{0}
h=ahess(g,0);
end
