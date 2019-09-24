%% Example 7: Curvature revisited

%%
% Examples 2 and 6 can be simplified using the programs |acurve| and |asurf|.

%%

% Example 2 revisited
t = ainit(linspace(0,2*pi,501),5);   
c = [cos(t).*cos(4*t); cos(t).*sin(4*t); cos(t)];
[kap,tau] = acurve(c);
figure(1), clf 
plot(t{0},kap{0},t{0},tau{0},t{0},kap{1},'--',t{0},tau{1},'--')
axis([0 2*pi -6 6]), grid on
title('Curvature, torsion, and their derivatives of a space curve')

% Example 6 revisited
[u,v] = ndgrid(linspace(-pi,pi));
[u,v] = ainit(u,v,2);
r = 2;
x = (r+cos(v/2)*sin(u)-sin(v/2)*sin(2*u))*cos(v);
y = (r+cos(v/2)*sin(u)-sin(v/2)*sin(2*u))*sin(v);
z = sin(v/2)*sin(u)+cos(v/2)*sin(2*u);
[K,H] = asurf([x;y;z]);
figure(2), clf, colormap(jet)
surf(x{0},y{0},z{0},K{0})
light, light, shading interp
axis equal, caxis([-1.5 .5]);
title('Klein bottle colored by Gaussian curvature')