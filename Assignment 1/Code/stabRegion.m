% Stability region computation example
clearvars
close all
clc

addpath(genpath('src\'));

% 1nd order Runge-Kutta growth factor
F1_2var = @(h,alpha) I + (h*A(alpha));
% 2nd order Runge-Kutta growth factor
F2_2var = @(h,alpha) I + (h*A(alpha)) + 0.5*(h*A(alpha)).^2;
% 4th order Runge-Kutta growth factor
F4_2var = @(h,alpha) I + (h*A(alpha)) + 1/2*(h*A(alpha)).^2 + ...
                1/6*(h*A(alpha)).^3 + 1/24*(h*A(alpha)).^4;

% alpha = deg2rad(0);
% h = 0.5;

I = diag(ones(2));
A = @(alpha) [0, 1; -1, 2.*cos(alpha)];
F1 = @(h) I + (h*A(alpha));
incr = 0.001;


alphaVec = linspace(0,pi,50);

%[solution,converge,info] = newton(f, xGuess, method ,toll, nmax , visualConfig)

for i=1:length(alphaVec)
    
    alpha = alphaVec(i);

    F1 = @(h) I + (h*A(alpha));
    prob = @(h) max(abs(eig( F1(h) ))) - 1;

    [h(i),conv(i)] = newton(prob, 2, 'f');
    
    
end

plot(alphaVec,h)

