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

alphalim = deg2rad([180 95]);

alphaVec = linspace(alphalim(1),alphalim(2),100);
guess = 1.1;
%[solution,converge,info] = newton(f, xGuess, method ,toll, nmax , visualConfig)

figure
hold on
axis equal
grid on

for i=1:length(alphaVec)
    
    alpha = alphaVec(i);
    %F4 = @(h) I + (h*A(alpha)) + 1/2*(h*A(alpha)).^2 + ...
    %            1/6*(h*A(alpha)).^3 + 1/24*(h*A(alpha)).^4;
    F2 = @(h) I + (h*A(alpha)) + 0.5*(h*A(alpha)).^2;
    F1 = @(h) I + (h*A(alpha));
    
    prob = @(h) max(abs(eig( F2(h) ))) - 1;
    
    %[h(i),conv(i)] = newton(prob, guess , 'f' ,1e-7, 1e4);
    [h(i),fval(i),conv(i)]=fsolve(prob,guess);
    
    eig_iterVec = eig(A(alpha));
    %eig_iterVec = eig(F1(h(i)));

    x(i) = h(i)*real(eig_iterVec(1));
    y(i) = h(i)*imag(eig_iterVec(1));
    
    guess = h(i);
    %if conv(i)==1
    plot(x(i),y(i),Marker="o",MarkerEdgeColor='k')
    % else
    % plot(x(i),y(i),Marker="o",MarkerEdgeColor='r')
    % end
end

%figure
%plot(cos(alphaVec),sin(alphaVec))

