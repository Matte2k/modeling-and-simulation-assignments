% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assugnment # 1
% Author: Matteo Baio 10667431

addpath(genpath('src\'));

%% Ex 1
%code ex1
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
axis padded
hold on

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


%%% RUNGE-KUTTA METHODS
% RK1
RK1 = rkSettings(1,[],[],[]);
[RK1.solution.x,RK1.solution.t,RK1.solution.info] = rungeKutta(f,x0,tmax,h,RK1);

% RK2
RK2 = rkSettings(2,'Heun',[],[]);
[RK2.solution.x,RK2.solution.t,RK2.solution.info] = rungeKutta(f,x0,tmax,h,RK2);

% RK4
RK4 = rkSettings(4,'Runge-Kutta',[],[]);
[RK4.solution.x,RK4.solution.t,RK4.solution.info] = rungeKutta(f,x0,tmax,h,RK4);


%%% ANALYTIC SOLUTION
t = RK1.solution.t;
fIVP_solution = @(t)(2.*t.^2 + 4.*t - exp(t) + 2);
plot(t, fIVP_solution(t), '--')
hold on

% %%% DEBUG
% fIVP_ode = @(t,x)(x - 2.*(t).^2 + 2);
% y0 = x_rk1(1);
% [tx, yx] = ode45(fIVP_ode, tspan, y0);
% plot(tx, yx, '--')

legend('rk1','rk2','rk4','ana','ode',Location='best')


%% EX BOH, debug startup






