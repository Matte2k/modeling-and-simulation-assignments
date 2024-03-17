% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assugnment # 1
% Author: Matteo Baio 10667431

addpath(genpath('src\'));

%% Ex 1
clearvars; close all; clc;

config = struct;
    config.print = true;
    config.plot = true;

f = @(x1,x2) [x2.^2 + x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
xGuess = [-1 1]';
toll = 1e-12;
nmax = 1e3;

figure()
grid on
hold on
axis padded

[solutionSym,convergeSym,infoSym] = newton(f, xGuess, 's', toll, nmax, config);
[solutionFD,convergeFD,infoFD] = newton(f, xGuess, 'f', toll, nmax, config);
[solutionCD,convergeCD,infoCD] = newton(f, xGuess, 'c', toll, nmax, config);
legend('sym','FD','CD',Location='best')

errSym = infoSym.errorVector(end);
errFD  = infoFD.errorVector(end);
errCD  = infoCD.errorVector(end);

if errSym < errFD && errSym < errCD
    fprintf('Sym is the best method\n');
elseif errFD < errSym && errFD < errCD
    fprintf('FD is the best method\n');
elseif errCD < errFD && errCD < errSym
    fprintf('CD is the best method\n');
end

%%% DEBUG - compare with built in function
% disp('----------------------------------------')
% fMatlab = @(x) [x(2).^2 + x(1) - 2; -x(1).^2 + x(2) + 10];
% fsolve(fMatlab, xGuess)


%% Ex 2
clearvars; close all; clc;

%%% DATA INPUT
tmax = 2; 
h = 0.5;                       % step size
x0 = 1;
f = @(x,t)(x - 2.*(t).^2 + 2);


%%% FIGURE INITIALIZATION
figure()
grid on
hold on
axis padded


%%% RK1
RK1 = struct;
    RK1.method = 'RK1';
    RK1.submethod = [];
    RK1.alpha = [];
    RK1.beta = [];
[x_rk1,t] = rungeKutta(f,x0,tmax,h,RK1);


%%% RK2
RK2 = struct;
    RK2.method = 'RK2';
    RK2.submethod = [];
    RK2.alpha = [];
    RK2.beta = [];
[x_rk2] = rungeKutta(f,x0,tmax,h,RK2);


%%% RK4
RK4 = struct;
    RK4.method = 'RK4';
    RK4.submethod = [];
    RK4.alpha = [];
    RK4.beta = [];
[x_rk4] = rungeKutta(f,x0,tmax,h,RK4);


%%% ANALYTIC SOLUTION
tspan = [0,tmax];
fIVP_solution = @(t)(2.*t.^2 + 4.*t - exp(t) + 2);
plot(t, fIVP_solution(t), '--')
hold on

% %%% DEBUG
% fIVP_ode = @(t,x)(x - 2.*(t).^2 + 2);
% y0 = x_rk1(1);
% [tx, yx] = ode45(fIVP_ode, tspan, y0);
% plot(tx, yx, '--')

legend('rk1','rk2','rk4','ana','ode',Location='best')
