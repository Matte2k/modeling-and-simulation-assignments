% Stability region computation example
clearvars
close all
clc

addpath(genpath('src\'));

% % 1nd order Runge-Kutta growth factor
% F1_2var = @(h,alpha) I + (h*A(alpha));
% % 2nd order Runge-Kutta growth factor
% F2_2var = @(h,alpha) I + (h*A(alpha)) + 0.5*(h*A(alpha)).^2;
% % 4th order Runge-Kutta growth factor
% F4_2var = @(h,alpha) I + (h*A(alpha)) + 1/2*(h*A(alpha)).^2 + ...
%                 1/6*(h*A(alpha)).^3 + 1/24*(h*A(alpha)).^4;

%%% Variables initialization
I = eye(2);
alphalim = deg2rad([179 91]);
alphaVec = linspace(alphalim(1),alphalim(2),50);
Afun = @(a) [0, 1; -1, 2*cos(a)];
guess = 2.4;  %rk1 = 2  rk2

% fsolve options settings
options = optimset('Display','iter');

% figure initialization
figure
hold on
axis equal
grid on

for i=1:length(alphaVec)
    
    alpha = alphaVec(i);
    A = Afun(alpha);

    %%% RK4 and RK2
    F4 = @(h) I + (h*A) + 0.5*(h*A).^2 + ...
                 1/6*(h*A).^3 + 1/24*(h*A).^4;  % - DON'T WORK -
    
    F2 = @(h) I + (h*A) + 0.5*(h*A).^2;         % - DON'T WORK -
    
    F1 = @(h) I + (h*A);
    
    prob = @(h) max(abs(eig(F1(h)))) - 1;    
    
    %abs(eig(F1(guess)))
    %fplot(prob,[0 2]);

    %%% H COMPUTATIONS
    %[hvec(i),conv(i)] = newton(prob, guess , 'f' ,1e-7, 1e3);
    %[hvec(i),fval(i),conv(i)]=fsolve(prob,guess);
    [hvec(i),fval(i),conv(i)]=fzero(prob,guess,options);
    
    guess = hvec(i);    % initial guess update

    %%% Eigenvalue of continuous problem
    eig_iterVec = eig(A);
    xA(i) = real(eig_iterVec(1));
    yA(i) = imag(eig_iterVec(1));
    plot(xA(i),yA(i),Marker="o",MarkerEdgeColor='k')
    
    %%% Eigenvalue of discrete problem
    xF(i) = hvec(i)*real(eig_iterVec(1));
    yF(i) = hvec(i)*imag(eig_iterVec(1));
    plot(xF(i),yF(i),Marker="+",MarkerEdgeColor='r')
    
end

% final stability region computed plot
figure(2)
plot(xF,yF)
axis equal
grid on


% KNOWN ISSUE:  with F2 and F4 the stability region is worng, the second
%               eigenvalue is way bigger than expected and so cause the h to be smaller
%               than expected.
