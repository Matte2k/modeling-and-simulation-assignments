% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assugnment # 1
% Author: Matteo Baio 10667431

addpath(genpath('src\'));

%% Ex 1
%code ex1
clearvars; close all; clc;

%%% DATA INPUT
f = @(x1,x2) [x2.^2 + x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
toll = 1e-12;
nmax = 1e3;
config = struct;
    config.print = true;
    config.plot = true;


%%% INITIAL GUESS DEFINITION
% Initial guess graphic analysis
f1 = @(x1,x2) x2.^2 + x1 - 2;
f2 = @(x1,x2) -x1.^2 + x2 + 10;
[Z1mesh,Z2mesh] = zerosGuess(f1,f2);
% -add graphic export here-
xGuess = [-2 -2]';

% Figure initialization
figure('Name','Convergence plot')
grid on
axis padded
hold on

%%% NEWTON'S METHODS
% Solution with different submethod
[solutionSym,convergeSym,infoSym] = newton(f, xGuess, 's', toll, nmax, config);
[solutionFD, convergeFD, infoFD]  = newton(f, xGuess, 'f', toll, nmax, config);
[solutionCD, convergeCD, infoCD]  = newton(f, xGuess, 'c', toll, nmax, config);
legend('sym','FD','CD',Location='best')
title('Newton convergence plot')
% -add graphic export here-

% Error compare between different submethod
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


%%% DEBUG
% compare with built in function
% disp('----------------------------------------')
% fMatlab = @(x) [x(2).^2 + x(1) - 2; -x(1).^2 + x(2) + 10];
% fsolve(fMatlab, xGuess)


%% Ex2 
%code ex2
clearvars; close all; clc;

%%% DATA INPUT
f = @(x,t)(x - 2.*(t).^2 + 2);
solutionIVP = @(t)(2.*t.^2 + 4.*t - exp(t) + 2);
x0   = 1;
tmax = 2; 
h = [0.5 0.2 0.05 0.01];     % step size

% Variable initialization
RK2 = cell(1,length(h));
RK4 = cell(1,length(h));
Analytic = cell(1,length(h));
lString  = cell(1,length(h)+1);


%%% ANALYTICAL SOLUTION & LEGEND STRING
for i = 1:length(h)
    Analytic{i}.t = 0:h(i):tmax;
    Analytic{i}.x = solutionIVP(Analytic{i}.t);
    lString{i} = ['h = ', num2str(h(i))];
end
lString{end} = 'Analytic';


%%% HEUN METHOD SOLUTION
% Figure initialization
figure('Name','RK2 results');
    t1 = tiledlayout(1,2);
    title(t1,'Heun method')

% Solution
nexttile
    hold on;    grid on;    box on;    axis padded
    for i = 1:length(h)                                         % RK2 solutions
        RK2{i} = rkSettings(2,'Heun');
        [RK2{i}.x,RK2{i}.t,RK2{i}.info] = rungeKutta(f,x0,tmax,h(i),RK2{i});
    end
    plot(Analytic{end}.t, Analytic{end}.x, '--','LineWidth',1); % Analytic solution
    title('Integration solution')
    legend(lString{1:end},'Location','best')
    xlabel('Time')
    ylabel('x')

% Error
nexttile       
    for i = 1:length(h)
        RK2{i}.error = Analytic{i}.x - RK2{i}.x;
        semilogy(RK2{i}.t,RK2{i}.error)
        hold on
    end
    grid on;    box on;    axis padded;
    title('Integration error')
    legend(lString{1:length(h)},'Location','best')
    xlabel('Time')
    ylabel('Error')

drawnow
% -add graphic export here-


%%% RK4 METHOD SOLUTION
% Figure initialization
figure('Name','RK4 results');
    t2 = tiledlayout(1,2);
    title(t2,'RK4 method')

% Solution
nexttile
    hold on;    grid on;    box on;    axis padded
    for i = 1:length(h)                                         % RK4 solutions
        RK4{i} = rkSettings(4);
        [RK4{i}.x,RK4{i}.t,RK4{i}.info] = rungeKutta(f,x0,tmax,h(i),RK4{i});
    end
    plot(Analytic{end}.t, Analytic{end}.x, '--','LineWidth',1); % Analytic solution
    title('Integration solution')
    legend(lString{1:end},'Location','best')
    xlabel('Time')
    ylabel('x')

% Error
nexttile       
    for i = 1:length(h)
        RK4{i}.error = Analytic{i}.x - RK4{i}.x;
        semilogy(RK4{i}.t,RK4{i}.error)
        hold on
    end
    grid on;    box on;    axis padded;
    title('Integration error')
    legend(lString{1:length(h)},'Location','best')
    xlabel('Time')
    ylabel('Error')
 
drawnow
% -add graphic export here-


%%% POST PROCESSING   //WIP
% Error vs CPU time
colorPalette = ['r','g','b','m'];           % Temporary
figure('Name','Trade off')
hold on;    grid on;    box on;     axis padded
for i = 1:length(h)
    scatter(RK2{i}.error(end),RK2{i}.info.timeCost,[],colorPalette(i),"filled","o")
    scatter(RK4{i}.error(end),RK4{i}.info.timeCost,[],colorPalette(i),"filled","square")
end
title('Integration solution')
legend(lString{1:end},'Location','best')    % Find better way to do
xlabel('Final error')
ylabel('Time cost')

drawnow
% -add graphic export here-


%%% DEBUG
% compare with built in function
% disp('----------------------------------------')
% fIVP_ode = @(t,x)(x - 2.*(t).^2 + 2);
% y0 = x0;
% [tx, yx] = ode45(fIVP_ode, [0 tmax], y0);
% plot(tx, yx, '--')
% 
% legend('rk1','rk2','rk4','ana','ode',Location='best')


%% Ex3
%code ex3
clearvars; close all; clc;


%% Ex4
%code ex4
clearvars; close all; clc;


%% Ex5
%code ex5
clearvars; close all; clc;


%% Ex6
%code ex6
clearvars; close all; clc;

%%% DATA INPUT
B = [-180.5, 219.5; 179.5, -220.5];             % IVP definition
f = @(x,t) [B(1,1).*x(1) + B(1,2).*x(2); ...
            B(2,1).*x(1) + B(2,2).*x(2)];       % function handle of RHS of the ODE
tmax = 1;       % start and ending time
x0 = [1 1]';    % initial guess
h = 0.1;        % time step

% Figure initialization
figure()
grid on
hold on
axis padded


%%% IVP SOLUTION
% RK4
RK4 = rkSettings(4,'Runge-Kutta');
[RK4.x,RK4.t,RK4.info] = rungeKutta(f,x0,tmax,h,RK4);

% IEX4
IEX4 = iSettings();
[IEX4.x,IEX4.t,IEX4.info] = iex4(f,x0,tmax,h,IEX4);

% Analytical    
solutionIVP =  @(t) expm(B.*t) * x0;            % 'expm' matrix exponential
Analytic.t = IEX4.t;
Analytic.x(:,1) = [1;1];
for i = 2 : (length(Analytic.t))       % CHIEDERE A RIC SE É IL MODO PIU' SMART
    tNow = Analytic.t(i);
    Analytic.x(:,i) = solutionIVP(tNow);
end
plot(Analytic.t,Analytic.x(1,:),'r')
plot(Analytic.t,Analytic.x(2,:),'r')


%%% EIGENVALUE ANALYSIS
%plot stability region for RK4 and IEX4
%plot eig(B)*h


%% Ex7
%code ex7
clearvars; close all; clc;

%%% DATA INPUT
f = @(x,t) [-5/2.*(1+8*sin(t)).*x(1); (1-x(1)).*x(2) + x(1)];
x0   = [1 1]';
tmax = 3;
h    = 0.1;

%%% AB
figure("Name",'Adams Bashforth')
hold on
grid on

abOptionsRK = abSettings(3,'RK');
[AB.startupRK.x,AB.startupRK.t,AB.startupRK.info] = ...
    adamsBashforth(f,x0,tmax,h,abOptionsRK);

abOptionsAB = abSettings(3,'AB');
[AB.startupAB.x,AB.startupAB.t,AB.startupAB.info] = ...
    adamsBashforth(f,x0,tmax,h,abOptionsAB);


%%% AM
figure("Name",'Adams Moulton')
hold on
grid on

amOptionsRK = amSettings(3,'RK');
[AM.startupRK.x,AM.startupRK.t,AM.startupRK.info] = ...
    adamsMoulton(f,x0,tmax,h,amOptionsRK);

amOptionsAM = amSettings(3,'AM');
[AM.startupAM.x,AM.startupAM.t,AM.startupAM.info] = ...
    adamsMoulton(f,x0,tmax,h,amOptionsAM);


%%% ABM
figure("Name",'Adams Bashforth Moulton')
hold on
grid on

abmOptionsRK = abmSettings(3,'RK');
[ABM.startupRK.x,ABM.startupRK.t,ABM.startupRK.info] = ...
    adamsBashforthMoulton(f,x0,tmax,h,abmOptionsRK);

abmOptionsABM = abmSettings(3,'ABM');
[ABM.startupABM.x,ABM.startupABM.t,ABM.startupABM.info] = ...
    adamsBashforthMoulton(f,x0,tmax,h,abmOptionsABM);


%%% BDF
figure("Name",'Backward Difference Formula')
hold on
grid on

bdfOptionsRK = bdfSettings(3,'RK');
[BDF.startupRK.x,BDF.startupRK.t,BDF.startupRK.info] =  ...
    backwardDifferenceFormula(f,x0,tmax,h,bdfOptionsRK);

bdfOptionsBDF = bdfSettings(3,'BDF');
[BDF.startupBDF.x,BDF.startupBDF.t,BDF.startupBDF.info] = ...
    backwardDifferenceFormula(f,x0,tmax,h,bdfOptionsBDF);


%%% COMPARE
figure("Name",'Multistep method compare')
hold on
grid on
plot(AB.startupRK.t,AB.startupRK.x,  'k-')
plot(AM.startupRK.t,AM.startupRK.x,  'r-')
plot(ABM.startupRK.t,ABM.startupRK.x,'g-')
plot(BDF.startupRK.t,BDF.startupRK.x,'b-')
xlabel('t')
ylabel('x')
legend('AM','AM','AB','AB','ABM','ABM','BDF','BDF','location','best')




