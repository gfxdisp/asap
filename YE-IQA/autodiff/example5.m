%% Example 5: Sinc function and rule of l'Hospital

%%

% Singularities of type "0/0" are resolved by l'Hospital's rule. 

% plot the sinc function and its derivatives
sinc = @(x) sin(x)./x;
x = ainit(linspace(-8,8),2);
y = sinc(x);
figure(1), clf
plot(x{0},y{0},x{0},y{1},x{0},y{2})
grid on, hold on

% mark values at multiples of pi
x = ainit(pi*(-2:2),2);
y = sinc(x);
plot(x{0},[y{0};y{1};y{2}],'ro')
[x{0};y{0};y{1};y{2}]
title('Sinc function and derivatives, accurately evaluated at x=0')
% Note that the second derivative at 0 is not plottet because 
% the order of x is 2 so that the second derivative cannot be
% determined there. Generally, when initializing x with order k,
% then derivatives up to order k-r are available, where r is
% the number of times you have to apply l'Hospital's rule until
% it yields a finite result.

% Accordingly, the order of the Taylor polynomial is reduced.
x = ainit(0,4);
disp('Quartic Taylor polynomial of sin(x):')
disp(ataylor(sin(x)))
disp('Cubic Taylor polynomial of sinc(x):')
disp(ataylor(sinc(x)))
disp('Quadratic Taylor polynomial of sinc(x^2):')
disp(ataylor(sinc(x.^2)))

% plot
t = linspace(-8,8);
figure(2), clf
plot(t,sinc(t),t,polyval(ataylor(sinc(x)),t))
grid on, axis([-8 8 -0.5 1])
title('Taylor polynomial of reduced degree')