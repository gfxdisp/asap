%% Example 2: Curvature and torsion of a space curve

%%

% initialize curve parameter
t = ainit(linspace(0,2*pi,501),3);               

% define curve as 3x1 audi array
c = [cos(t)*cos(4*t); cos(t)*sin(4*t); cos(t)];

% compute curvature kap and torsion tau 
b   = cross(c{1},c{2});
kap = sqrt(sum(b.^2)./sum(c{1}.^2).^3);   
tau = dot(b,c{3})./sum(b.^2); 

% plot curve and a few tangents 
figure(1), clf     
D = aget(c,[1:10:51 226:10:276]);             % select data
T = [D{0};D{0}+.2*D{1}];                      % tangents
plot3(c{0}(1,:),c{0}(2,:),c{0}(3,:)), hold on
plot3(T([1 4],:),T([2 5],:),T([3 6],:),'r')
grid on, axis equal
title('Space curve and a few tangents')

% plot
figure(2), clf
plot(t{0},[kap;tau])
grid on, axis tight
title('Curvature and torsion')