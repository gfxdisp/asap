%% Example 4: Laplacian of the peaks function and its maximum

%%

% initialize audi grid and evaluate
[u,v] = ndgrid(linspace(-3,3,50));
[u,v] = ainit(u,v,2);                        
P = peaks(u,v);
L = alap(P);

% plot peaks function and its Laplacian
figure(1), clf
mesh(u{0},v{0},P{0})
axis tight, grid on  
title('Peaks function')
figure(2), clf
mesh(u{0},v{0},L{0})
axis tight, grid on, hold on
title('Laplacian and its maximum')

% Search the maximum of the Laplacian using Newton's method 
% for the gradient of L starting form [u,v] = [0.0,-1.5].
% This requests the Hessian of L, and hence fourth derivatives
% of the given function.
[u,v] = ainit(0.0,-1.5,4);
w = [u;v];
for i = 1:4
  L = alap(peaks(w(1),w(2)));
  w = w - ahess(L,0)\agrad(L,0);
end

% At the resulting point u=w(1),v=w(2), the norm of the 
% gradient of the Laplacian is less than 1e-12.
plot3(w{0}(1),w{0}(2),L{0},'r*')