%% Example 6: Gaussian and mean curvature of a Klein bottle

%%

% initialize audi grid
[u,v] = ndgrid(linspace(-pi,pi));
[u,v] = ainit(u,v,2);

% define parametrization of a Klein bottle
r = 2;
x = (r+cos(v/2)*sin(u)-sin(v/2)*sin(2*u))*cos(v);
y = (r+cos(v/2)*sin(u)-sin(v/2)*sin(2*u))*sin(v);
z = sin(v/2)*sin(u)+cos(v/2)*sin(2*u);
S = [x;y;z];

% differential geometry
J = ajac(S);                              % Jacobian
N = cross(J(:,1),J(:,2));                 % normal vector
N = N/norm(N);                            % normalize
G = J'*J;                                 % first fundamental form
B = [N'*adiff(S,2,0) N'*adiff(S,1,1);...  % second fundamental form
     N'*adiff(S,1,1) N'*adiff(S,0,2)];
W = G\B;                                  % Weingarten map
K = det(W);                               % Gaussian curvature
H = trace(W)/2;                           % mean curvature

% color surface by Gaussian curvature
figure(1), clf, colormap(jet)
surf(x{0},y{0},z{0},K{0})
light, light, shading interp
axis equal, caxis([-1.5 .5]);
title('Gaussian curvature')

% color surface by mean curvature
figure(2), clf, colormap(jet)
surf(x{0},y{0},z{0},H{0})
light, light, shading interp
axis equal, caxis([-1 1])
title('Mean curvature, discontinuous since surface is not orientable')