% This function uses automatic differentiation to compute the gradient of 
% an anonymous function f(x) at the point x.
function grad=quick_agrad(f,x)
a=num2cell(x);
% Must assign ainit variables at same time.
adims=size(a);
N=adims(2);
a2=cell(1,N);
[a2{:}]=ainit(a{:},2);

% Define auto-diff version of f(x) for the value given in a.
g=f([a2{:}]);

% Evaluates the Hessian at the values given in a (the 0 means that we
% take 0th order derivatives of the gradient (i.e. the gradient itself)).
% This is equivalent to:
% grad=agrad(g)
% grad{0}
grad=agrad(g,0);
end
