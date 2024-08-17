% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment # 1
% Author: Matteo Baio 232805

% Written and executed on MATLAB 2023a - Windows 10
addpath(genpath('src/'));
addpath(genpath('figure/'));


%% Ex 1
%code ex1
clearvars; close all; clc;
graphicSettings;

%%% DATA INPUT
f = @(x1,x2) [x2.^2 - x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
xGuess = [[3 2]',[3 -2]'];
toll = 1e-6;
nmax = 1e3;
avgTimeIter = 50;
config = struct;
    config.print = true;
    config.plot = false;

%%% INITIAL GUESS DEFINITION
% Initial guess graphic analysis
f1 = @(x1,x2) x2.^2 - x1 - 2;           f2 = @(x1,x2) -x1.^2 + x2 + 10;
[Z1mesh,Z2mesh] = zerosGuess(f1,f2);            % plot functions
plot(xGuess(1,:),xGuess(2,:),'ko');             % plot guess
legend('$f_1=0$','$f_2=0$','Initial guess')     % plot legend
set(gcf,'units','centimeters','position',[0,0,20,9]);
exportgraphics(gcf,'figure\ex1_initGuess.eps');
%exportgraphics(gcf,'figure\ex1_initGuess.png',Resolution=1500);

% Variable initialization
numZeros = size(xGuess,2);
solutionSym = zeros(2,numZeros);    solutionFD  = zeros(2,numZeros);    solutionCD  = zeros(2,numZeros);
convergeSym = zeros(1,numZeros);    convergeFD =  zeros(1,numZeros);    convergeCD =  zeros(1,numZeros);
infoSym = cell(1,numZeros);         infoFD = cell(1,numZeros);          infoCD = cell(1,numZeros);

%%% NEWTON'S METHODS
for i = 1:numZeros
    [solutionSym(:,i),convergeSym(i),infoSym{i}] = newton(f, xGuess(:,i), 's', toll, nmax, avgTimeIter, config); % symbolic math
    [solutionFD(:,i), convergeFD(i), infoFD{i}]  = newton(f, xGuess(:,i), 'f', toll, nmax, avgTimeIter, config); % forward difference
    [solutionCD(:,i), convergeCD(i), infoCD{i}]  = newton(f, xGuess(:,i), 'c', toll, nmax, avgTimeIter, config); % centered difference

    % Error compare between different submethod
    errSym(i) = infoSym{i}.absError;
    errFD(i)  = infoFD{i}.absError;
    errCD(i)  = infoCD{i}.absError;
    
end

%%% RESULTS PLOT
figure('Name','Convergence zero'); 
t1 = tiledlayout(1,2);
     tileTitle = title(t1,'Newton''s error evolution');     tileTitle.Interpreter = 'latex';
for i=1:numZeros
    nexttile
        semilogy((1:infoFD{i}.iteration),infoFD{i}.errorVector,'o-');   % forward diff error plot
        hold on;    grid on;    axis padded;    box on; 
        semilogy((1:infoCD{i}.iteration),infoCD{i}.errorVector,'--');   % centered diff error plot
        semilogy((1:infoSym{i}.iteration),infoSym{i}.errorVector,'*');  % symolic math error plot
        yl = yline(toll,'--','Label','tollerance');                             % tollerance value
            yl.LabelHorizontalAlignment = 'left';
            yl.Interpreter = 'latex';
        xlabel('$Iteration$');    ylabel('$Error$');
        ylim([10e-12,10e0])
        legend('$FD$','$CD$','$sym$',Location='northeast')
        tString = sprintf('$z_{%d}=[%.4f,%.4f]$',i,solutionFD(1,i),solutionFD(2,i));
        title(tString)
end
set(gcf,'units','centimeters','position',[0,0,20,9]);
exportgraphics(gcf,'figure\ex1_convergence.eps');


%% Ex2 
%code ex2
clearvars; close all; clc;
[cmap]=graphicSettings();

%%% DATA INPUT
f = @(x,t)(x - 2.*(t).^2 + 2);
solutionIVP = @(t)(2.*t.^2 + 4.*t - exp(t) + 2);
avgTimeIter = 1000;
x0   = 1;
tmax = 2; 
h = [0.5 0.2 0.05 0.01]';     % step size

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
    tileTitle = title(t1,'Heun method');
    tileTitle.Interpreter = 'latex';

% Solution
nexttile
    hold on;    grid on;    box on;    axis padded
    for i = 1:length(h)                                         % RK2 solutions
        RK2{i} = rkSettings(2,'Heun',[],[],avgTimeIter);
        [RK2{i}.x,RK2{i}.t,RK2{i}.info] = rungeKutta(f,x0,tmax,h(i),RK2{i});
    end
    plot(Analytic{end}.t, Analytic{end}.x, '--','LineWidth',1); % Analytic solution
    title('Integration solution')
    legend(lString{1:end},'Location','best')
    xlabel('$Time$');   ylabel('$x$');

% Error
nexttile       
    for i = 1:length(h)
        RK2{i}.error = Analytic{i}.x - RK2{i}.x;
        semilogy(RK2{i}.t,RK2{i}.error,'o-',MarkerSize=5)
        hold on
    end
    grid on;    box on;    axis padded;
    title('Integration error')
    legend(lString{1:length(h)},'Location','best')
    xlabel('$Time$');   ylabel('$Error$');

drawnow
set(gcf,'units','centimeters','position',[0,0,20,11]);
exportgraphics(gcf,'figure\ex2_heun.eps');

%%% RK4 METHOD SOLUTION
% Figure initialization
figure('Name','RK4 results');
    t2 = tiledlayout(1,2);
    tileTitle = title(t2,'RK4 method');
    tileTitle.Interpreter = 'latex';

% Solution
nexttile
    hold on;    grid on;    box on;    axis padded
    for i = 1:length(h)                                         % RK4 solutions
        RK4{i} = rkSettings(4,'Runge-Kutta',[],[],avgTimeIter);
        [RK4{i}.x,RK4{i}.t,RK4{i}.info] = rungeKutta(f,x0,tmax,h(i),RK4{i});
    end
    plot(Analytic{end}.t, Analytic{end}.x, '--','LineWidth',1); % Analytic solution
    title('Integration solution')
    legend(lString{1:end},'Location','best')
    xlabel('$Time$');   ylabel('$x$');

% Error
nexttile       
    for i = 1:length(h)
        RK4{i}.error = Analytic{i}.x - RK4{i}.x;
        semilogy(RK4{i}.t,RK4{i}.error,'o-',MarkerSize=5)
        hold on
    end
    grid on;    box on;    axis padded;
    title('Integration error')
    legend(lString{1:length(h)},'Location','best')
    xlabel('$Time$');   ylabel('$Error$'); 
drawnow
set(gcf,'units','centimeters','position',[0,0,20,11]);
exportgraphics(gcf,'figure\ex2_rk4.eps');

%%% POST PROCESSING
% Error vs CPU time
figure('Name','Trade off')
cmap = cmap(1:length(h),:);             % color map definition
    cbTicksCount = 1:length(h);
    cbTicksPos = [0.5, cbTicksCount, (cbTicksCount(end)+0.5)];
hold on;    grid on;    box on;     axis padded
for i = 1:length(h)
    scatter(RK2{i}.error(end),RK2{i}.info.avgTimeCost,50,cmap(i,:),"filled","o");
    scatter(RK4{i}.error(end),RK4{i}.info.avgTimeCost,60,cmap(i,:),"filled","square")
end

colormap(cmap)                              % apply colormap
clim([cbTicksPos(1),cbTicksPos(end)])       
    cb = colorbar;                          % apply colorbar
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$h$';
    cb.Ticks = cbTicksPos;
    cb.TickLabels = {'',num2str(h),''};
set(cb,'TickLabelInterpreter','latex')

title('Error vs. Cost comparison');      legend('RK2','RK4','Location','best');     % title + legend
xlabel('Final error');   ylabel('Time cost');                                       % axis label

drawnow
set(gcf,'units','centimeters','position',[0,0,12,10]);
exportgraphics(gcf,'figure\ex2_error.eps');


%% Ex3
%code ex3
clearvars; close all; clc;
graphicSettings;

dimSys = 2;
stepNum = 100;

%%% RK2 STABILITY PROBLEM
figure('Name','RK2');
    t1 = tiledlayout(1,2);
    tileTitle = title(t1,'RK2 stability analysis');
    tileTitle.Interpreter = 'latex';

nexttile
    [RK2.F,RK2.guess,RK2.degVec,RK2.hVec]=selectOp('RK2',dimSys);
    stabGuess(RK2.F,RK2.degVec(1),RK2.hVec);

nexttile
    [RK2.xF,RK2.yF,RK2.hvec,RK2.stabEdge]=stabRegion(RK2.F,stepNum,RK2.degVec,RK2.guess);
    xlim([-4 4]);   ylim([-4 4]);
    
leg = legend('$RK2$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
set(gcf,'units','centimeters','position',[0,0,20,9]);
exportgraphics(gcf,'figure\ex3_rk2.eps');

fprintf('Solution of the problem in alpha = pi using RK2 is:\n h = %.4f \n\n',RK2.hvec(1));

%%% RK4 STABILITY PROBLEM
figure('Name','RK4');
    t2 = tiledlayout(1,2);
    tileTitle = title(t2,'RK4 stability analysis');
    tileTitle.Interpreter = 'latex';

nexttile
    [RK4.F,RK4.guess,RK4.degVec,RK4.hVec]=selectOp('RK4',dimSys);
    stabGuess(RK4.F,RK4.degVec(1),RK4.hVec);
    
nexttile
    [RK4.xF,RK4.yF,RK4.hvec,RK4.stabEdge]=stabRegion(RK4.F,stepNum,RK4.degVec,RK4.guess);
    xlim([-4 4]);   ylim([-4 4]);

leg = legend('$RK4$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
set(gcf,'units','centimeters','position',[0,0,20,9]);
exportgraphics(gcf,'figure\ex3_rk4.eps');

fprintf('Solution of the problem in alpha = pi using RK4 is:\n h = %.4f \n\n',RK4.hvec(1));

%%% STABILITY COMPARE RK2 vs RK4
figure('Name','RK4 vs. RK2');
    axis equal;     grid on;    box on;     hold on
    stabRK4 = fill(RK4.stabEdge(1,:),RK4.stabEdge(2,:),'r','FaceAlpha',0.3);
        stabRK4.LineWidth = 1.5;
        stabRK4.EdgeColor = 'r';
    stabRK2 = fill(RK2.stabEdge(1,:),RK2.stabEdge(2,:),'b','FaceAlpha',0.3);
        stabRK2.LineWidth = 1.5;
        stabRK2.EdgeColor = 'b';
    xline(0,LineStyle='--');        yline(0,LineStyle='--');
    xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
    ylim([-3.5 3.5]);
    legend('$RK4$','$RK2$')
    title('$RK2$ vs $RK4$')
set(gcf,'units','centimeters','position',[0,0,12,10]);
exportgraphics(gcf,'figure\ex3_stability.eps');
%exportgraphics(gcf,'figure\ex3_stability.png',Resolution=1500);


%% Ex4
%code ex4
clearvars; close all; clc;
[cmap]=graphicSettings();

%%% INPUT
x0      = [1 1]';
tmax    = 1;
tollVec = [1e-3 1e-4 1e-5 1e-6]';
degVec  = [180 0];
dimSys  = 2;
stepNum = 100;
orders  = [1 2 4];

RKcell  = cell(1,length(orders));
[RKcell{1}]=selectOp('RK1',dimSys);
[RKcell{2}]=selectOp('RK2',dimSys);
[RKcell{3}]=selectOp('RK4',dimSys);

Afun = @(a) [0, 1; -1, 2*cos(a)];
alphalim = deg2rad(degVec);
alphaVec = linspace(alphalim(1),alphalim(2),stepNum);

%%% VARIABLES INITIALIZATION
Gcell    = cell(length(orders),length(tollVec));
solEdge = cell(length(orders),length(tollVec));
hVec    = zeros(length(orders),length(tollVec));
hStart  = zeros(length(orders),length(tollVec));
feval = zeros(length(orders),length(tollVec));
xF    = zeros(1,length(alphaVec));
yF    = zeros(1,length(alphaVec));

%%% COMPUTE SOLUTION RK1,2,4 
for ord = 1:length(RKcell) 
    RK = RKcell{ord};
    method = sprintf('RK%.0f',orders(ord));
    
    for tols = 1:length(tollVec)
        guess = [tollVec(tols) 1];          % Initial guess definition
    
        for i=1:length(alphaVec)
            alpha = alphaVec(i);
            A = Afun(alpha);
    
            %%% PROBLEM DEFINITION
            nstep = @(h) (tmax-0)/h;
            rkSol = @(h) (RK(h,A))^nstep(h) * x0;   % Runge-Kutta solution
            analSol = expm(A*tmax)*x0;              % Analytic solution
            
            G = @(h) norm((analSol-rkSol(h)),inf)-tollVec(tols);    % Problem function   
            [hVec(ord,tols),~,conv]=fzero(G,guess);                 % Problem solution

            if i == 1
                Gcell{ord,tols} = G;                                   % Initial guess analysis
                feval(ord,tols) = orders(ord) * nstep(hVec(ord,tols)); % Function evaluation count
                hStart(ord,tols) = hVec(ord,tols);
                fprintf('Solution of the problem in alpha = pi using %s and toll = %.4e is:\n h = %.4e \n\n',method,tollVec(tols),hStart(ord,tols));
            end
            guess = hVec(ord,tols);    % initial guess update
                   
            % Eigenvalue of discrete problem
            eig_iterVec = eig(A);
            xF(i) = hVec(ord,tols)*real(eig_iterVec(1));
            yF(i) = hVec(ord,tols)*imag(eig_iterVec(1));

            if conv < 1
                warning('Implicit equation at the step %d has not been solved correctly',i)
            end
        end
        
        solEdge{ord,tols} = [xF,flip(xF,2);yF,-flip(yF,2)];    % Solution edge
    end
end

%%% PLOT SOLUTION RK1,2,4
% Visual settings
hMax = 0.5;
view = struct;
    view.ax1 = struct;
    view.ax2 = struct;

view.ax1.x = [-1e-3 2.5e-3  ; -0.01   0.06   ; -0.05   0.45  ];
view.ax1.y = [-2e-3   2e-3  ; -1.5e-4 0.5e-4 ; -1.5e-3 0.5e-3];
view.ax2.x = [-2.5e-3 3e-3  ; -0.07   0.11   ; -0.6 0.9];
view.ax2.y = [-2e-3   2e-3  ; -0.09   0.09   ; -0.7 0.7];

% Plot loop
for ord = 1:length(RKcell)
    tstring = sprintf('RK%.0f',orders(ord));
    figure("Name",tstring)
    t = tiledlayout(1,2);
    tileTitle = title(t,tstring);
    tileTitle.Interpreter = 'latex';
    
    ax1 = nexttile;
        hold (ax1,'on')
        grid on;    axis padded;    box on;
        yline(0, LineStyle = '--')
        xlim(view.ax1.x(ord,:));    ylim(view.ax1.y(ord,:));
        xlabel('$h$');              ylabel('$G(h)$');
        title('Stability function in $\alpha = \pi$')

    ax2 = nexttile;
        hold (ax2,'on')
        grid on;    axis padded;    box on;
        xline(0,'--');              yline(0,'--')
        xlim(view.ax2.x(ord,:));    ylim(view.ax2.y(ord,:));
        xlabel('$Re\{h\lambda\}$'); ylabel('$Im\{h\lambda\}$');
        title('Stability region')
    
    cmap = cmap(1:length(tollVec),:);             % color map definition
        cbTicksCount = 1:length(tollVec);
        cbTicksPos = [0.5, cbTicksCount, (cbTicksCount(end)+0.5)];
    
    for tols = 1:length(tollVec)
        fplot(ax1,Gcell{ord,tols},[0 hMax],'Color',cmap(tols,:),LineWidth=1.2);
        xline(ax1,hVec(ord,tols),'Color',cmap(tols,:),'LineStyle',':');
        plot (ax2,solEdge{ord,tols}(1,:), solEdge{ord,tols}(2,:));
    end
    
    ax1.Layer = 'top';      hold(ax1,'off');
    ax2.Layer = 'top';      hold(ax2,'off');
    
    colormap(cmap)                              % apply colormap
    clim([cbTicksPos(1),cbTicksPos(end)])       
        cb = colorbar;                          % apply colorbar
        cb.Label.Interpreter = 'latex';
        cb.Label.String = '$toll$';
        cb.Ticks = cbTicksPos;
        cb.TickLabels = {'',num2str(tollVec),''};
    set(cb,'TickLabelInterpreter','latex')
    
    figName = sprintf('figure/ex4_%s.eps',tstring);
    set(gcf,'units','centimeters','position',[0,0,21,9]);
    exportgraphics(gcf,figName);
end

%%% ZOOM RK1 PLOT
% Visual settings
view.ax1.xZoom = [-1e-5  3e-5];     view.ax1.yZoom = [-6e-5  6e-5];
view.ax2.xZoom = [-2.5e-5 3e-5];    view.ax2.yZoom = [-2e-5   2e-5];
ord = 1;

% Plot loop
tstring = sprintf('RK1 zoom');
figure("Name",tstring)
t = tiledlayout(1,2);
tileTitle = title(t,tstring);
tileTitle.Interpreter = 'latex';
    
    ax1 = nexttile;
        hold (ax1,'on')
        grid on;    axis padded;    box on;
        yline(0, LineStyle = '--')
        xlim(view.ax1.xZoom(ord,:));    ylim(view.ax1.yZoom(ord,:));
        xlabel('$h$');                  ylabel('$G(h)$');
        title('Stability function in $\alpha = \pi$')

    ax2 = nexttile;
        hold (ax2,'on')
        grid on;    axis padded;    box on;
        xline(0,'--');                  yline(0,'--')
        xlim(view.ax2.xZoom(ord,:));    ylim(view.ax2.yZoom(ord,:));
        xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
        title ('Stability region')
    
    cmap = cmap(1:length(tollVec),:);             % color map definition
        cbTicksCount = 1:length(tollVec);
        cbTicksPos = [0.5, cbTicksCount, (cbTicksCount(end)+0.5)];
    
    for tols = 1:length(tollVec)
        fplot(ax1,Gcell{ord,tols},[0 hMax],'Color',cmap(tols,:),LineWidth=1.2);
        xline(ax1,hVec(ord,tols),'Color',cmap(tols,:),'LineStyle',':');
        plot (ax2,solEdge{ord,tols}(1,:), solEdge{ord,tols}(2,:));
    end

    ax1.Layer = 'top';      hold(ax1,'off');
    ax2.Layer = 'top';      hold(ax2,'off');

    colormap(cmap)                              % apply colormap
    clim([cbTicksPos(1),cbTicksPos(end)])       
        cb = colorbar;                          % apply colorbar
        cb.Label.Interpreter = 'latex';
        cb.Label.String = '$toll$';
        cb.Ticks = cbTicksPos;
        cb.TickLabels = {'',num2str(tollVec),''};
    set(cb,'TickLabelInterpreter','latex')
    set(gcf,'units','centimeters','position',[0,0,21,9]);
    exportgraphics(gcf,'figure/ex4_RK1zoom.eps');

%%% FUNCTION EVALUATION vs TOLLERANCE
figure("Name",'Feval at alpha=pi')
lString  = cell(1,length(RKcell));
for ord = 1:length(RKcell) 
    loglog(tollVec,feval(ord,:),'o-')
    lString{ord} = ['RK', num2str(orders(ord))];
    hold on
end
grid on;    axis padded;    box on;
xlabel('Tol');      ylabel('Function evaluations')
title('Function evaluations for $\alpha = \pi$')
legend(lString)
set(gcf,'units','centimeters','position',[0,0,9,7]);
exportgraphics(gcf,'figure/ex4_feval.eps');


%% Ex5
%code ex5
clearvars; close all; clc;
[cmap]=graphicSettings();

thetaVec = [0.1 0.3 0.4 0.7 0.9]';
dimSys   = 2;
stepNum  = 100;
solEdge = cell(1,length(thetaVec));
lstring  = cell(1,length(thetaVec));
pal = 'grmcb';

%%% BI2 STABILITY PROBLEM
for i = 1:length(thetaVec)
    figure('Name','BI2');
        t = tiledlayout(1,2);
        tstring = sprintf('$BI2_{%.1f}$ stability analysis',thetaVec(i));
        tileTitle = title(t,tstring);
        tileTitle.Interpreter = 'latex';
    
    nexttile
        [F,guess,degVec,hVec]=selectOp('BI2',dimSys,thetaVec(i));
        stabGuess(F,degVec(1),hVec);
    
    nexttile
        [~,~,~,solEdge{i}]=stabRegion(F,stepNum,degVec,guess);
        xlim([-5.5 10.5]);
    
    lstring{i} = sprintf('$BI2_{%.1f}$',thetaVec(i));
    leg = legend(lstring{i},'$Stable$','Orientation', 'Horizontal');
        leg.Layout.Tile = 'south';        leg.Color = 'w';
    figName = sprintf('figure/ex5_BI2_%.1f.eps',thetaVec(i));
    set(gcf,'units','centimeters','position',[0,0,11,9]);
    exportgraphics(gcf,figName);
end

figure('Name','BI2 compare');
    axis equal;     grid on;    box on;     hold on
    cmap = cmap(1:length(thetaVec),:);             % color map definition
        cbTicksCount = 1:length(thetaVec);
        cbTicksPos = [0.5, cbTicksCount, (cbTicksCount(end)+0.5)];
        
    for i = 1:length(thetaVec)
        plot(solEdge{i}(1,:),solEdge{i}(2,:),'Color',cmap(i,:))
    end

    colormap(cmap)                              % apply colormap
    clim([cbTicksPos(1),cbTicksPos(end)])       
        cb = colorbar;                          % apply colorbar
        cb.Label.Interpreter = 'latex';
        cb.Label.String = '$theta$';
        cb.Ticks = cbTicksPos;
        cb.TickLabels = {'',num2str(thetaVec),''};
    set(cb,'TickLabelInterpreter','latex')

xline(0,LineStyle='--');        yline(0,LineStyle='--');
xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
xlim([-6 11]);
tstring = sprintf('$BI2_{n}$ stability analysis');
title(tstring)
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex5_stabReg.eps');
    

%% Ex6
%code ex6
clearvars; close all; clc;
[cmap]=graphicSettings();

%%% DATA INPUT
B = [-180.5, 219.5; 179.5, -220.5];             % IVP definition
f = @(x,t) [B(1,1).*x(1) + B(1,2).*x(2); ...
            B(2,1).*x(1) + B(2,2).*x(2)];       % function handle of RHS of the ODE
tmax = 5;       % start and ending time
x0 = [1 1]';    % initial guess
h  = 0.1;       % time step

%%% IVP SOLUTION
graphicOutput = false;
% RK4
RK4 = rkSettings(4,'Runge-Kutta');
[RK4.x,RK4.t,RK4.info] = rungeKutta(f,x0,tmax,h,RK4,graphicOutput);

% RK4mod
RK4mod = rkSettings(4,'Runge-Kutta');   hmod = 0.005;
[RK4mod.x,RK4mod.t,RK4mod.info] = rungeKutta(f,x0,tmax,hmod,RK4,graphicOutput);

% IEX4
IEX4 = iSettings();
[IEX4.x,IEX4.t,IEX4.info] = iex4(f,x0,tmax,h,IEX4,graphicOutput);

% Analytical    
solutionIVP =  @(t) expm(B.*t) * x0;            % 'expm' matrix exponential
Analytic.t = IEX4.t;
Analytic.x(:,1) = [1;1];
for i = 2 : (length(Analytic.t))
    tNow = Analytic.t(i);
    Analytic.x(:,i) = solutionIVP(tNow);
end

%%% INTEGRATION RESULTS
figure('Name','Integration results');
    t1 = tiledlayout(1,2);
    tileTitle = title(t1,'Integration results');
    tileTitle.Interpreter = 'latex';

nexttile
    axis padded;     grid on;    box on;     hold on
    plot(RK4.t,RK4.x(1,:),'Color',cmap(1,:))
    plot(IEX4.t,IEX4.x(1,:),'Color',cmap(2,:))
    %plot(RK4mod.t,RK4mod.x(1,:),'--','Color',cmap(3,:))
    plot(Analytic.t,Analytic.x(1,:),'o','Color',cmap(5,:),MarkerSize=2)
    xlabel('$Time$');   ylabel('$x_1$');
    title('$x_1$ integration')
    ylim([-0.2 1.2])
    
nexttile
    axis padded;     grid on;    box on;     hold on
    plot(RK4.t,RK4.x(2,:),'Color',cmap(1,:))
    plot(IEX4.t,IEX4.x(2,:),'Color',cmap(2,:))
    %plot(RK4mod.t,RK4mod.x(2,:),'--','Color',cmap(3,:))
    plot(Analytic.t,Analytic.x(2,:),'o','Color',cmap(5,:),MarkerSize=2)
    xlabel('$Time$');   ylabel('$x_2$');
    title('$x_2$ integration')
    ylim([-0.2 1.2])

leg = legend('$RK4$','$IEX4$','$Analytical$','Orientation', 'Horizontal');
%leg = legend('$RK4$','$IEX4$','$RK4mod$','$Analytical$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
set(gcf,'units','centimeters','position',[0,0,20,8]);
exportgraphics(gcf,'figure/ex6_integOutput.eps');    
%exportgraphics(gcf,'figure/ex6_integOutputMod.eps');

%%% EIGENVALUE ANALYSIS
dimSys  = length(x0);
stepNum = 100;
eigB = eig(B);
    xEigs = h*real(eigB);           yEigs = h*imag(eigB);
    xEigsMod = hmod*real(eigB);     yEigsMod = hmod*imag(eigB);   

% IEX4
figure('Name','IEX4');
    t1 = tiledlayout(1,2);
    tileTitle = title(t1,'IEX4 stability analysis');
    tileTitle.Interpreter = 'latex';

nexttile
    [IEX4.F,IEX4.stability.guess,IEX4.stability.degVec,IEX4.stability.hMax]=selectOp('IEX4',dimSys);
    stabGuess(IEX4.F,IEX4.stability.degVec(1),IEX4.stability.hMax);

nexttile
    [IEX4.stability.xF,IEX4.stability.yF,IEX4.stability.hvec,IEX4.stability.stabEdge]=stabRegion(IEX4.F,stepNum,IEX4.stability.degVec,IEX4.stability.guess);
    xlim([-5 15])
    
leg = legend('$IEX4$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';      leg.Color = 'w';
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex6_iex4.eps');

% RK4    
figure('Name','RK4');
    t2 = tiledlayout(1,2);
    tileTitle = title(t2,'RK4 stability analysis');
    tileTitle.Interpreter = 'latex';

nexttile
    [RK4.F,RK4.guess,RK4.degVec,RK4.hVec]=selectOp('RK4',dimSys);
    stabGuess(RK4.F,RK4.degVec(1),RK4.hVec);

nexttile
    [RK4.xF,RK4.yF,RK4.hvec,RK4.stabEdge]=stabRegion(RK4.F,stepNum,RK4.degVec,RK4.guess);
    xlim([-4 1])

leg = legend('$RK4$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex6_rk4.eps');

% IEX4 vs RK4 vs matrixB
figure('Name','IEX4 vs. RK2');
    axis equal;     grid on;    box on;     hold on;
    stabIEX4 = fill(IEX4.stability.stabEdge(1,:),IEX4.stability.stabEdge(2,:),'b','FaceAlpha',0.3);
        stabIEX4.LineWidth = 1.5;
        stabIEX4.EdgeColor = 'b';

    stabRK4 = fill(RK4.stabEdge(1,:),RK4.stabEdge(2,:),'r','FaceAlpha',0.3);
        stabRK4.LineWidth = 1.5;
        stabRK4.EdgeColor = 'r';
    
    eigsPlot = plot(xEigs,yEigs);
        eigsPlot.LineStyle = "none";
        eigsPlot.Marker = "o";
        eigsPlot.MarkerEdgeColor = 'k';
        eigsPlot.MarkerFaceColor = [0.50 0.50 0.50];
        eigsPlot.MarkerSize = 5;
        
    eigsPlotMod = plot(xEigsMod,yEigsMod);
        eigsPlotMod.LineStyle = "none";
        eigsPlotMod.Marker = "square";
        eigsPlotMod.MarkerEdgeColor = 'k';
        eigsPlotMod.MarkerFaceColor = [0.50 0.50 0.50];
        eigsPlotMod.MarkerSize = 5;
    
    xline(0,LineStyle='--');       yline(0,LineStyle='--');
    xlabel('$Re_{h\lambda}$');     ylabel('$Im_{h\lambda}$');
    xlim([-45 15]);                ylim([-10 10]);
    legend('$IEX4$','$RK4$','$h\lambda$','$h_{mod}\lambda$',Orientation='Horizontal',Location='southoutside')
    title('$RK4$ vs $IEX4$')
set(gcf,'units','centimeters','position',[0,0,11,6]);
exportgraphics(gcf,'figure/ex6_stabReg.eps');
%exportgraphics(gcf,'figure/ex6_stabReg.png',Resolution=1500);


%% Ex7
%code ex7
clearvars; close all; clc;
[cmap]=graphicSettings();

%%% DATA INPUT
f = @(x,t) [-5/2*(1+8*sin(t))*x(1); (1-x(1))*x(2) + x(1)];
eig1 = @(t) -5/2*(1+8*sin(t));
eig2 = @(t) 1-exp(-5/2*t+20*(cos(t)-1));
x0     = [1 1]';
tmax   = 3;
h      = 0.1;
visual = false;

%%% INTEGRATION
abOptionsRK = abSettings(3,'RK');
[AB.startupRK.x,AB.startupRK.t,AB.startupRK.info] = ...
    adamsBashforth(f,x0,tmax,h,abOptionsRK,visual);

abOptionsAB = abSettings(3,'AB');
[AB.startupAB.x,AB.startupAB.t,AB.startupAB.info] = ...
    adamsBashforth(f,x0,tmax,h,abOptionsAB,visual);

amOptionsRK = amSettings(3,'RK');
[AM.startupRK.x,AM.startupRK.t,AM.startupRK.info] = ...
    adamsMoulton(f,x0,tmax,h,amOptionsRK,visual);

amOptionsAM = amSettings(3,'AM');
[AM.startupAM.x,AM.startupAM.t,AM.startupAM.info] = ...
    adamsMoulton(f,x0,tmax,h,amOptionsAM,visual);

abmOptionsRK = abmSettings(3,'RK');
[ABM.startupRK.x,ABM.startupRK.t,ABM.startupRK.info] = ...
    adamsBashforthMoulton(f,x0,tmax,h,abmOptionsRK,visual);

abmOptionsABM = abmSettings(3,'ABM');
[ABM.startupABM.x,ABM.startupABM.t,ABM.startupABM.info] = ...
    adamsBashforthMoulton(f,x0,tmax,h,abmOptionsABM,visual);

bdfOptionsRK = bdfSettings(3,'RK');
[BDF.startupRK.x,BDF.startupRK.t,BDF.startupRK.info] =  ...
    backwardDifferenceFormula(f,x0,tmax,h,bdfOptionsRK,visual);

bdfOptionsBDF = bdfSettings(3,'BDF');
[BDF.startupBDF.x,BDF.startupBDF.t,BDF.startupBDF.info] = ...
    backwardDifferenceFormula(f,x0,tmax,h,bdfOptionsBDF,visual);


%%% PLOT
figure("Name",'Adams Bashforth')
t1 = tiledlayout(2,1);
    tileTitle = title(t1,'Adams Bashforth');
    tileTitle.Interpreter = 'latex';

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(AB.startupRK.t, AB.startupRK.x(1,:), 'b');
    plot(AB.startupAB.t, AB.startupAB.x(1,:), 'r--');
    ylabel('$x_{1}$');      xlabel('Time');
    ylim([-10 10]);

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(AB.startupRK.t, AB.startupRK.x(2,:), 'b');
    plot(AB.startupAB.t, AB.startupAB.x(2,:), 'r--');
    ylabel('$x_{2}$');      xlabel('Time');
    ylim([-10 10]);
    legend('Runge-Kutta starter', 'Adams Bashforth starter', ...
                Location='southoutside', Orientation='horizontal');
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex7_ab.eps');

% AM
figure("Name",'Adams Moulton')
t1 = tiledlayout(2,1);
    tileTitle = title(t1,'Adams Moulton');
    tileTitle.Interpreter = 'latex';

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(AM.startupRK.t, AM.startupRK.x(1,:), 'b');
    plot(AM.startupAM.t, AM.startupAM.x(1,:), 'r--');
    ylabel('$x_{1}$');      xlabel('Time');

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(AM.startupRK.t, AM.startupRK.x(2,:), 'b');
    plot(AM.startupAM.t, AM.startupAM.x(2,:), 'r--');
    ylabel('$x_{2}$');      xlabel('Time');
    legend('Runge-Kutta starter', 'Adams Moulton starter', ...
                Location='southoutside', Orientation='horizontal');
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex7_am.eps');

% ABM
figure("Name",'Adams Bashforth Moulton')
t1 = tiledlayout(2,1);
    tileTitle = title(t1,'Adams Bashforth Moulton');
    tileTitle.Interpreter = 'latex';

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ABM.startupRK.t,  ABM.startupRK.x(1,:),  'b');
    plot(ABM.startupABM.t, ABM.startupABM.x(1,:), 'r--');
    ylabel('$x_{1}$');      xlabel('Time');

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ABM.startupRK.t,  ABM.startupRK.x(2,:),  'b');
    plot(ABM.startupABM.t, ABM.startupABM.x(2,:), 'r--');
    ylabel('$x_{2}$');      xlabel('Time');
    legend('Runge-Kutta starter', 'Adams Bashforth Moulton starter', ...
                Location='southoutside', Orientation='horizontal');
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex7_abm.eps');

% BDF
figure("Name",'Backward Difference Formula')
t1 = tiledlayout(2,1);
    tileTitle = title(t1,'Backward Difference Formula');
    tileTitle.Interpreter = 'latex';

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(BDF.startupRK.t,  BDF.startupRK.x(1,:),  'b');
    plot(BDF.startupBDF.t, BDF.startupBDF.x(1,:), 'r--');
    ylabel('$x_{1}$');      xlabel('Time');

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(BDF.startupRK.t,  BDF.startupRK.x(2,:),  'b');
    plot(BDF.startupBDF.t, BDF.startupBDF.x(2,:), 'r--');
    ylabel('$x_{2}$');      xlabel('Time');
    legend('Runge-Kutta starter', 'Backward Difference Formula starter', ...
                Location='southoutside', Orientation='horizontal');
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex7_bdf.eps');

%%% COMPARE
figure("Name",'Multistep method compare')
t1 = tiledlayout(2,1);
    tileTitle = title(t1,'Multistep method compare');
    tileTitle.Interpreter = 'latex';

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(AB.startupRK.t,  AB.startupRK.x(1,:),  'Color', cmap(1,:))
    plot(AM.startupRK.t,  AM.startupRK.x(1,:),  'Color', cmap(2,:))
    plot(ABM.startupRK.t, ABM.startupRK.x(1,:), 'Color', cmap(3,:))
    plot(BDF.startupRK.t, BDF.startupRK.x(1,:), 'o', 'Color', cmap(5,:), MarkerSize=3)
    ylabel('$x_{1}$');      xlabel('Time');
    ylim([-4 2]);

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(AB.startupRK.t,  AB.startupRK.x(2,:),  'Color', cmap(1,:))
    plot(AM.startupRK.t,  AM.startupRK.x(2,:),  'Color', cmap(2,:))
    plot(ABM.startupRK.t, ABM.startupRK.x(2,:), 'Color', cmap(3,:))
    plot(BDF.startupRK.t, BDF.startupRK.x(2,:), 'o', 'Color', cmap(5,:), MarkerSize=3)
    ylabel('$x_{2}$');      xlabel('Time');
    ylim([-2 22]);
    legend('AM','AB','ABM','BDF', Location='southoutside', Orientation='horizontal');
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex7_compare.eps');


%%% EIGENVALUE ENVELOP
timeVec = 0:h:tmax;
figure("Name",'Eigenvalue evolution')
    hold on;    grid on;    axis padded;    box on;
    plot(timeVec,h*eig1(timeVec),'r')
    plot(timeVec,h*eig2(timeVec),'b')
    limitAB3 = yline(-0.6,Color=cmap(1,:),LineWidth=1.5,LineStyle='--',Label='AB3');
        limitAB3.LabelHorizontalAlignment='center';
        limitAB3.Interpreter='latex';
    limitAM3 = yline(-6,Color=cmap(2,:),LineWidth=1.5,LineStyle='--',Label='AM3');
        limitAM3.LabelHorizontalAlignment='center';
        limitAM3.Interpreter='latex';
    limitABM3 = yline(-1.7,Color=cmap(5,:),LineWidth=1.5,LineStyle='--',Label='ABM3');
        limitABM3.LabelHorizontalAlignment='center';
        limitABM3.Interpreter='latex';
    xlabel('Time');     ylabel('$h\lambda_{i}$');     ylim([-6.75,0.75])
    legend('$h\lambda_{1}$','$h\lambda_{2}$',Location='best')
set(gcf,'units','centimeters','position',[0,0,11,9]);
exportgraphics(gcf,'figure/ex7_eig.eps');

%%% --- END CODE ---


%% FUNCTIONS

function [x,t,info] = ab1(f,x0,tmax,h,abOptions)
%AB1 - Adams Bashforth method of order 1
%
%   Syntax:
%       [x,t,info] = ab1(f,x0,tmax,h,abOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abOtpions(*),   struct:  see abSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abOptions.method = [];
        abOptions.alpha  = [];
        abOptions.beta   = [];
    end

    %%% Default method
    if isempty(abOptions.method)
        abOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abOptions.method
        case 'Standard'
            alpha = 1;      % standard alpha parameter
            beta  = 1;      % standard beta parameter

            if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))                                % TO DEBUG
                warning('Parameters matrix unused, %s parameters are used instead',...
                    abOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abOptions.alpha) || isempty(abOptions.beta)                                          % TO DEBUG
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(abOptions.alpha),[1,1])) || not(isequal(size(abOptions.beta),[1,1]))    % TO DEBUG
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = abOptions.alpha;    % custom alpha parameter
            beta = abOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x =  [x0 zeros(dimSys,length(t)-1)];          % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% AB1 loop
    for i = 1 : (length(t)-1)       % main loop of the method
        fk = fvalVec(:,i);
        x(:,i+1) = x(:,i) + h/alpha * (beta*fk);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;     % function evaluation counter update
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = ab2(f,x0,tmax,h,abOptions)
%AB2 - Adams Bashforth method of order 2
%
%   Syntax:
%       [x,t,info] = ab2(f,x0,tmax,h,abOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,2]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abOtpions(*),   struct:  see abSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abOptions.method = [];
        abOptions.alpha  = [];
        abOptions.beta   = [];
    end

    %%% Default method
    if isempty(abOptions.method)
        abOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abOptions.method
        case 'Standard'
            alpha = 2;          % standard alpha parameter
            beta = [3 -1];      % standard beta parameter

            if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))
                warning('Parameters matrix unused, %s parameters are used instead',...
                    abOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abOptions.alpha) || isempty(abOptions.beta)
                error('No parameters matrix has been given as input');    % missing parameters matrix

            elseif not(isequal(size(abOptions.alpha),[1,1])) || not(isequal(size(abOptions.beta),[1,2]))
                error('Parameters matrix dimensions are invalid');        % parameters matrix with wrong size

            end

            alpha = abOptions.alpha;    % custom alpha parameter
            beta = abOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-2)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), ...
                    zeros(dimSys,length(t)-2)];     % fval vector allocation
    feval = feval + 2*dimSys;                       % function evaluation counter update

    %%% AB2 loop
    for i = 2 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        x(:,i+1) = x(:,i) + h/alpha * (beta(1)*fk1 + beta(2)*fk2);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;     % function evaluation counter update
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = ab3(f,x0,tmax,h,abOptions)
%AB3 - Adams Bashforth method of order 3
%
%   Syntax:
%       [x,t,info] = ab3(f,x0,tmax,h,abOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,3]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abOtpions(*),   struct:  see abSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abOptions.method = [];
        abOptions.alpha  = [];
        abOptions.beta   = [];
    end

    %%% Default method
    if isempty(abOptions.method)
        abOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abOptions.method
        case 'Standard'
            alpha = 12;             % standard alpha parameter
            beta = [23 -16  5];     % standard beta parameter

            if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))
                warning('Parameters matrix unused, %s parameters are used instead',...
                    abOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abOptions.alpha) || isempty(abOptions.beta)
                error('No parameters matrix has been given as input');    % missing parameters matrix

            elseif not(isequal(size(abOptions.alpha),[1,1])) || not(isequal(size(abOptions.beta),[1,3]))
                error('Parameters matrix dimensions are invalid');        % parameters matrix with wrong size

            end

            alpha = abOptions.alpha;    % custom alpha parameter
            beta = abOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-3)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), ...
                    zeros(dimSys,length(t)-3)];     % fval vector allocation
    feval = feval + 3*dimSys;                       % function evaluation counter update

    %%% AB3 loop
    for i = 3 : (length(t)-1)        % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        x(:,i+1) = x(:,i) + h/alpha * (beta(1)*fk1 + beta(2)*fk2 + beta(3)*fk3);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;      % function evaluation counter updateend
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;               % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = ab4(f,x0,tmax,h,abOptions)
%AB4 - Adams Bashforth method of order 4
%
%   Syntax:
%       [x,t,info] = ab4(f,x0,tmax,h,abOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,4]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abOtpions(*),   struct:  see abSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abOptions.method = [];
        abOptions.alpha  = [];
        abOptions.beta   = [];
    end

    %%% Default method
    if isempty(abOptions.method)
        abOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abOptions.method
        case 'Standard'
            alpha = 24;                 % standard alpha parameter
            beta = [55 -59  37 -9];     % standard beta parameter

            if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))
                warning('Parameters matrix unused, %s parameters are used instead',...
                    abOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abOptions.alpha) || isempty(abOptions.beta)
                error('No parameters matrix has been given as input');    % missing parameters matrix

            elseif not(isequal(size(abOptions.alpha),[1,1])) || not(isequal(size(abOptions.beta),[1,4]))
                error('Parameters matrix dimensions are invalid');        % parameters matrix with wrong size

            end

            alpha = abOptions.alpha;    % custom alpha parameter
            beta = abOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-4)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), f(x(:,4),t(4)), ...
                    zeros(dimSys,length(t)-4)];     % fval vector allocation
    feval = feval + 4*dimSys;                       % function evaluation counter update

    %%% AB4 loop
    for i = 4 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        fk4 = fvalVec(:,i-3);
        x(:,i+1) = x(:,i) + h/alpha * (beta(1)*fk1 + beta(2)*fk2 + beta(3)*fk3 + beta(4)*fk4);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;     % function evaluation counter updateend
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = abm1(f,x0,tmax,h,abmOptions)
%ABM1 - Adams Bashforth Moulton method of order 1
%
%   Syntax:
%       [x,t,info] = abm1(f,x0,tmax,h,abmOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions(*),  struct:  see abmSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abmOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abmOptions.method = [];
        abmOptions.alpha  = [];
        abmOptions.betaP  = [];
        abmOptions.betaC  = [];
    end

    %%% Default method
    if isempty(abmOptions.method)
        abmOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abmOptions.method
        case 'Standard'
            alpha = 1;      % standard alpha parameter
            betaP = 1;      % standard beta predictor parameter
            betaC = 1;      % standard beta corrector parameter

            if not(isempty(abmOptions.alpha)) || not(isempty(abmOptions.betaP)) || not(isempty(abmOptions.betaC))G
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    abmOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abmOptions.alpha) || isempty(abmOptions.betaP) || isempty(abmOptions.betaC)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(abmOptions.alpha),[1,1])) || not(isequal(size(abmOptions.betaP),[1,1])) || not(isequal(size(abmOptions.betaC),[1,1]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = abmOptions.alpha;     % custom alpha parameter
            betaP = abmOptions.betaP;     % custom beta predictor parameter
            betaC = abmOptions.betaC;     % beta corrector parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)]; % fval vector allocation
    feval = feval + dimSys;                     % function evaluation counter update

    %%% ABM1 loop
    for i = 1 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        xp = x(:,i) + h/alpha(1) * (betaP(1)*fk1);  
        tp = t(i) + h;
        x(:,i+1) = x(:,i) + h/alpha(1) * (betaC(1)*f(xp,tp));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;   % function evaluation counter update
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = abm2(f,x0,tmax,h,abmOptions)
%ABM2 - Adams Bashforth Moulton method of order 2
%
%   Syntax:
%       [x,t,info] = abm2(f,x0,tmax,h,abmOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,2]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions(*),  struct:  see abmSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abmOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abmOptions.method = [];
        abmOptions.alpha  = [];
        abmOptions.betaP  = [];
        abmOptions.betaC  = [];
    end

    %%% Default method
    if isempty(abmOptions.method)
        abmOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abmOptions.method
        case 'Standard'
            alpha = 2;          % standard alpha parameter
            betaP = [3 -1];     % standard beta predictor parameter
            betaC = [1  1];     % standard beta corrector parameter

            if not(isempty(abmOptions.alpha)) || not(isempty(abmOptions.betaP)) || not(isempty(abmOptions.betaC))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    abmOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abmOptions.alpha) || isempty(abmOptions.betaP) || isempty(abmOptions.betaC)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(abmOptions.alpha),[1,1])) || not(isequal(size(abmOptions.betaP),[1,2])) || not(isequal(size(abmOptions.betaC),[1,2]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = abmOptions.alpha;     % custom alpha parameter
            betaP = abmOptions.betaP;     % custom beta predictor parameter
            betaC = abmOptions.betaC;     % beta corrector parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-2)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), ...
                    zeros(dimSys,length(t)-2)]; % fval vector allocation
    feval = feval + 2*dimSys;                   % function evaluation counter update

    %%% ABM2 loop
    for i = 2 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        xp = x(:,i) + h/alpha * (betaP(1)*fk1 + betaP(2)*fk2);
        tp = t(i) + h;
        x(:,i+1) = x(:,i) + h/alpha * (betaC(1)*f(xp,tp) + betaC(2)*fk1);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;   % function evaluation counter update
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = abm3(f,x0,tmax,h,abmOptions)
%ABM3 - Adams Bashforth Moulton method of order 3
%
%   Syntax:
%       [x,t,info] = abm3(f,x0,tmax,h,abmOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,3]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions(*),  struct:  see abmSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abmOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abmOptions.method = [];
        abmOptions.alpha  = [];
        abmOptions.betaP  = [];
        abmOptions.betaC  = [];
    end

    %%% Default method
    if isempty(abmOptions.method)
        abmOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abmOptions.method
        case 'Standard'
            alpha = 12;             % standard alpha parameter
            betaP = [23 -16  5];    % standard beta predictor parameter
            betaC = [5   8  -1];    % standard beta corrector parameter

            if not(isempty(abmOptions.alpha)) || not(isempty(abmOptions.betaP)) || not(isempty(abmOptions.betaC))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    abmOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abmOptions.alpha) || isempty(abmOptions.betaP) || isempty(abmOptions.betaC)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(abmOptions.alpha),[1,1])) || not(isequal(size(abmOptions.betaP),[1,3])) || not(isequal(size(abmOptions.betaC),[1,3]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = abmOptions.alpha;     % custom alpha parameter
            betaP = abmOptions.betaP;     % custom beta predictor parameter
            betaC = abmOptions.betaC;     % beta corrector parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-3)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), ...
                    zeros(dimSys,length(t)-3)];     % fval vector allocation
    feval = feval + 3*dimSys;                       % function evaluation counter update

    %%% ABM3 loop
    for i = 3 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        xp = x(:,i) + h/alpha * (betaP(1)*fk1 + betaP(2)*fk2 + betaP(3)*fk3);
        tp = t(i) + h;
        x(:,i+1) = x(:,i) + h/alpha * (betaC(1)*f(xp,tp) + betaC(2)*fk1 + betaC(3)*fk2);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;   % function evaluation counter updateend
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = abm4(f,x0,tmax,h,abmOptions)
%ABM4 - Adams Bashforth Moulton method of order 4
%
%   Syntax:
%       [x,t,info] = abm3(f,x0,tmax,h,abmOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,4]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions(*),  struct:  see abmSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       abmOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abmOptions.method = [];
        abmOptions.alpha  = [];
        abmOptions.betaP  = [];
        abmOptions.betaC  = [];
    end

    %%% Default method
    if isempty(abmOptions.method)
        abmOptions.method = 'Standard';       % Heun parameters sets as default
    end

    %%% Parameters definition
    switch abmOptions.method
        case 'Standard'
            alpha = 24;                 % standard alpha parameter
            betaP = [55 -59  37  -9];   % standard beta predictor parameter
            betaC = [9   19  -5   1];   % standard beta corrector parameter

            if not(isempty(abmOptions.alpha)) || not(isempty(abmOptions.betaP)) || not(isempty(abmOptions.betaC))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    abmOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abmOptions.alpha) || isempty(abmOptions.betaP) || isempty(abmOptions.betaC)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(abmOptions.alpha),[1,1])) || not(isequal(size(abmOptions.betaP),[1,4])) || not(isequal(size(abmOptions.betaC),[1,4]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = abmOptions.alpha;     % custom alpha parameter
            betaP = abmOptions.betaP;     % custom beta predictor parameter
            betaC = abmOptions.betaC;     % beta corrector parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-4)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), f(x(:,4),t(4)), ...
                    zeros(dimSys,length(t)-4)];     % fval vector allocation
    feval = feval + 4*dimSys;                       % function evaluation counter update

    %%% ABM4 loop
    for i = 4 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        fk4 = fvalVec(:,i-3);
        xp = x(:,i) + h/alpha * (betaP(1)*fk1 + betaP(2)*fk2 + betaP(3)*fk3 + betaP(4)*fk4);
        tp = t(i) + h;
        x(:,i+1) = x(:,i) + h/alpha * (betaC(1)*f(xp,tp) + betaC(2)*fk1 + betaC(3)*fk2 + betaC(4)*fk3);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;   % function evaluation counter updateend
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [abmOptions] = abmSettings(order,startup,method,alpha,betaP,betaC)
%AB MSETTINGS - Create the options struct for Adams Bashforth Moulton methods
%
%   Syntax:
%       [abmOptions] = abmSettings(order,startup,method,alpha,betaP,betaC)
%
%   Input:
%       order,          double:  select ABM method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom' 
%       alpha(*),  double[1,1]:  insert custom parameter alpha for ABM method
%       betaP(*),  double[1,n]:  insert custom parameter beta predictor for ABM method
%       betaP(*),  double[1,n]:  insert custom parameter beta corrector for ABM method
%
%   Output:
%       abmOptions,  struct:  contains settings for AB method:
%           - abmOptions.order = order
%           - abmOptions.startup = startup;
%           - abmOptions.method = method
%           - abmOptions.alpha = alpha
%           - abmOptions.betaP = betaP
%           - abmOptions.betaC = betaC
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside abm# function
%       alpha:   set as empty by default
%       betaP:    set as empty by default
%       betaC:    set as empty by default
%


    %%% Default value for optional input
    if nargin < 6    
        betaC = [];
        if nargin < 5
            betaP = [];
            if nargin < 4          
                alpha = [];
                if nargin < 3
                    method = [];
                    if nargin < 2
                        startup = [];
                    end
                end
            end
        end
    end

    %%% Struct definition
    abmOptions = struct;
        abmOptions.order   = order;
        abmOptions.startup = startup;
        abmOptions.method  = method;
        abmOptions.alpha   = alpha;
        abmOptions.betaP   = betaP;
        abmOptions.betaC   = betaC;

end

function [abOptions] = abSettings(order,startup,method,alpha,beta)
%AB SETTINGS - Create the options struct for Adams Bashforth methods
%
%   Syntax:
%       [abOptions] = abSettings(order,startup,method,alpha,beta)
%
%   Input:
%       order,          double:  select AB method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom' 
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,n]:  insert custom parameter beta for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AB method:
%           - abOptions.order = order
%           - abOptions.startup = startup;
%           - abOptions.method = method
%           - abOptions.alpha = alpha
%           - abOptions.beta = beta
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside ab# function
%       alpha:   set as empty by default
%       beta:    set as empty by default
%


    %%% Default value for optional input
    if nargin < 5
        beta = [];
        if nargin < 4
            alpha = [];
            if nargin < 3 
                method = [];
                if nargin < 2
                    startup = [];
                end
            end
        end
    end

    %%% Struct definition
    abOptions = struct;
        abOptions.order   = order;
        abOptions.method  = method;
        abOptions.alpha   = alpha;
        abOptions.beta    = beta;
        abOptions.startup = startup;

end

function [x,t,info] = adamsBashforth(f,x0,tmax,h,abOptions,visualConfig)
%ADAMS BASHFORTH - Adams Bashforth method selection and application
%
%   Syntax:
%       [x,t,info] = adamsBashforth(f,x0,tmax,h,abOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abOtpions,      struct:  see abSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    %%% Dimension Check for initial guess
    if size(x0,2) > 1 && visualConfig == true
        % warning to remember that initial guess correct input format
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step',h);
    end

    if size(x0,2) > abOptions.order
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected')

    elseif size(x0,2) < abOptions.order 
        if isempty(abOptions.startup)
            error('Need to define startup method to retrive a valid initial guess')
        end
    else
        infoStartup = [];   % initialization of the infoStartup variable used in info
        if not(isempty(abOptions.startup))
            warning('Startup method defined will not be used, x0 is already valid')
        end
    end

    %%% Adams Bashforth order selector based on 'abOptions'
    switch abOptions.order
        case 1
            [x,t,infoMethod] = ab1(f,x0,tmax,h,abOptions);

        case 2
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 2
                guessOrder = 2;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abOptions);
            end
            
            [x,t,infoMethod] = ab2(f,x0,tmax,h,abOptions);
            
        case 3
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 3
                guessOrder = 3;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abOptions);
            end

            [x,t,infoMethod] = ab3(f,x0,tmax,h,abOptions);

        case 4
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 4
                guessOrder = 4;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abOptions);
            end

            [x,t,infoMethod] = ab4(f,x0,tmax,h,abOptions);

        otherwise
            error('Please insert a valid order as input');
    end

    if visualConfig == true
        plot(t,x,'-');     % plot of the solution
    end

    if nargout == 3
        info = struct;      % info struct build-up
            info.fvalVec   = infoMethod.fvalVec;        

        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost  = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
            
            % update the implict flag checking the startup method implicity
            if infoStartup.implicit == true
                info.implicit = true;
            else
                info.implicit = false;
            end
        else
            info.timeCost  = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
            info.implicit  = false;
        end
    end

end

function [x,t,info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptions,visualConfig)
%ADAMS BASHFORTH MOULTON - Adams Bashforth Moulton method selection and application
%
%   Syntax:
%       [x,t,info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions,      struct:  see abmSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    %%% Dimension Check for initial guess
    if size(x0,2) > 1 && visualConfig == true
        % warning to remember that initial guess correct input format
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step',h);
    end

    if size(x0,2) > abmOptions.order
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')

    elseif size(x0,2) < abmOptions.order 
        if isempty(abmOptions.startup)
            error('Need to define startup method to retrive a valid initial guess')
        end

    else
        infoStartup = [];   % initialization of the infoStartup variable used in info
        if not(isempty(abmOptions.startup))
            warning('Startup method defined will not be used, x0 is already valid')
        end
    end

    %%% Adams Bashforth Moulton order selector based on 'abmOptions'
    switch abmOptions.order
        case 1
            [x,t,infoMethod] = abm1(f,x0,tmax,h,abmOptions);

        case 2
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 2
                guessOrder = 2;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abmOptions);
            end
            
            [x,t,infoMethod] = abm2(f,x0,tmax,h,abmOptions);
            
        case 3
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 3
                guessOrder = 3;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abmOptions);
            end

            [x,t,infoMethod] = abm3(f,x0,tmax,h,abmOptions);

        case 4
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 4
                guessOrder = 4;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abmOptions);
            end

            [x,t,infoMethod] = abm4(f,x0,tmax,h,abmOptions);

        otherwise
            error('Please insert a valid order as input');
    end

    if visualConfig == true
        plot(t,x,'-');     % plot of the solution
    end

    if nargout == 3
        info = struct;      % info struct build-up
            info.fvalVec   = infoMethod.fvalVec;        

        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost  = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
            
            % update the implict flag checking the startup method implicity
            if infoStartup.implicit == true
                info.implicit = true;
            else
                info.implicit = false;
            end
        else
            info.timeCost  = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
            info.implicit  = false;
        end
    end

end

function [x,t,info] = adamsMoulton(f,x0,tmax,h,amOptions,visualConfig)
%ADAMS MOULTON - Adams Moulton method selection and application
%
%   Syntax:
%       [x,t,info] = adamsMoulton(f,x0,tmax,h,amOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions,      struct:  see amSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end
               
    %%% Dimension Check for initial guess
    if size(x0,2) > 1 && visualConfig == true
        % warning to remember that initial guess correct input format
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step',h);
    end

    % case for AM order grater than 1
    if amOptions.order > 1
        if size(x0,2) > amOptions.order-1 
            error('The initial guess is invalid, too many input in x0 vector compare to the order selected')

        elseif size(x0,2) < amOptions.order-1 
            if isempty(amOptions.startup)
                error('Need to define startup method to retrive a valid initial guess')
            end

        else
            infoStartup = [];   % initialization of the infoStartup variable used in info
            if not(isempty(amOptions.startup))
                warning('Startup method defined will not be used, x0 is already valid')
            end
        end
    end

    % case for AM order equal to 1
    if amOptions.order == 1
        if size(x0,2) > amOptions.order 
            error('The initial guess is invalid, too many input in x0 vector compare to the order selected')
        
        elseif size(x0,2) < amOptions.order 
            if isempty(amOptions.startup)
                error('Need to define startup method to retrive a valid initial guess')
            end
        
        else
            infoStartup = [];   % initialization of the infoStartup variable used in info
            if not(isempty(amOptions.startup))
                warning('Startup method defined will not be used, x0 is already valid')
            end
        end
    end

    %%% Adams Bashforth order selector based on 'amOptions'
    switch amOptions.order
        case 1
            [x,t,infoMethod] = am1(f,x0,tmax,h,amOptions);

        case 2
            [x,t,infoMethod] = am2(f,x0,tmax,h,amOptions);
            
        case 3
            % retrive a valid initial guess based on 'amOptions'
            if size(x0,2) < 2
                guessOrder = 2;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,amOptions);
            end

            [x,t,infoMethod] = am3(f,x0,tmax,h,amOptions);

        case 4
            % retrive a valid initial guess based on 'amOptions'
            if size(x0,2) < 3
                guessOrder = 3;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,amOptions);
            end

            [x,t,infoMethod] = am4(f,x0,tmax,h,amOptions);

        otherwise
            error('Please insert a valid order as input');
    end

    if visualConfig == true
        plot(t,x,'-');     % plot of the solution
    end

    if nargout == 3
        info = struct;
            info.fvalVec  = infoMethod.fvalVec; 
            info.implicit = true;

        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost  = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
        else
            info.timeCost  = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
        end
    end

end

function [x,t,info] = am1(f,x0,tmax,h,amOptions)
%AM1 - Adams Moulton method of order 1
%
%   Syntax:
%       [x,t,info] = am1(f,x0,tmax,h,amOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions(*),   struct:  see amSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       amOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        amOptions.method  = [];
        amOptions.options = [];
        amOptions.alpha   = [];
        amOptions.beta    = [];
    end

    %%% Default method
    if isempty(amOptions.method)
        amOptions.method = 'Standard';
    end

    if isempty(amOptions.options)
        amOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );    % default fsolve options
    end

    %%% Parameters definition
    switch amOptions.method
        case 'Standard'
            alpha = 1;      % standard alpha parameter
            beta = 1;       % standard beta parameter

            if not(isempty(amOptions.alpha)) || not(isempty(amOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    amOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(amOptions.alpha) || isempty(amOptions.beta)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(amOptions.alpha),[1,1])) || not(isequal(size(amOptions.beta),[1,1]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = amOptions.alpha;    % custom alpha parameter
            beta = amOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)]; % fval vector allocation
    feval = feval + dimSys;                     % function evaluation counter update

    %%% AM1 loop
    for i = 1 : (length(t)-1)         % main loop of the method
        tp = t(i) + h;
        fp = @(xp) x(:,i) + h/alpha(1) * (beta(1)*f(xp,tp)) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), amOptions.options);
        feval = feval + info.funcCount;

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);    % timer stop

    if nargout == 3
        info = struct;                % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [x,t,info] = am2(f,x0,tmax,h,amOptions)
%AM2 - Adams Moulton method of order 2
%
%   Syntax:
%       [x,t,info] = am2(f,x0,tmax,h,amOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions(*),   struct:  see amSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       amOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        amOptions.method  = [];
        amOptions.options = [];
        amOptions.alpha   = [];
        amOptions.beta    = [];
    end

    %%% Default method
    if isempty(amOptions.method)
        amOptions.method = 'Standard';
    end

    if isempty(amOptions.options)
        amOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );    % default fsolve options
    end

    %%% Parameters definition
    switch amOptions.method
        case 'Standard'
            alpha = 2;      % standard alpha parameter
            beta = [1 1];   % standard beta parameter

            if not(isempty(amOptions.alpha)) || not(isempty(amOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    amOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(amOptions.alpha) || isempty(amOptions.beta)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(amOptions.alpha),[1,1])) || not(isequal(size(amOptions.beta),[1,2]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = amOptions.alpha;    % custom alpha parameter
            beta = amOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];           % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% AM2 loop
    for i = 1 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        tp = t(i) + h;
        fp = @(xp) x(:,i) + h/alpha * (beta(1)*f(xp,tp) + beta(2)*fk1) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), amOptions.options);
        feval = feval + info.funcCount;

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [x,t,info] = am3(f,x0,tmax,h,amOptions)
%AM3 - Adams Moulton method of order 3
%
%   Syntax:
%       [x,t,info] = am3(f,x0,tmax,h,amOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,2]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions(*),   struct:  see amSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       amOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        amOptions.method  = [];
        amOptions.options = [];
        amOptions.alpha   = [];
        amOptions.beta    = [];
    end

    %%% Default method
    if isempty(amOptions.method)
        amOptions.method = 'Standard';
    end

    if isempty(amOptions.options)
        amOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );    % default fsolve options
    end

    %%% Parameters definition
    switch amOptions.method
        case 'Standard'
            alpha = 12;         % standard alpha parameter
            beta = [5 8 -1];    % standard beta parameter

            if not(isempty(amOptions.alpha)) || not(isempty(amOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    amOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(amOptions.alpha) || isempty(amOptions.beta)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(amOptions.alpha),[1,1])) || not(isequal(size(amOptions.beta),[1,3]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = amOptions.alpha;    % custom alpha parameter
            beta = amOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-2)];           % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), ...
                    zeros(dimSys,length(t)-2)];   % fval vector allocation
    feval = feval + 2*dimSys;                     % function evaluation counter update

    %%% AM3 loop
    for i = 2 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        tp = t(i) + h;
        fp = @(xp) x(:,i) + h/alpha * (beta(1)*f(xp,tp) + beta(2)*fk1 + beta(3)*fk2) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), amOptions.options);
        feval = feval + info.funcCount; 

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [x,t,info] = am4(f,x0,tmax,h,amOptions)
%AM4 - Adams Moulton method of order 4
%
%   Syntax:
%       [x,t,info] = am4(f,x0,tmax,h,amOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,3]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions(*),   struct:  see amSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       amOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        amOptions.method  = [];
        amOptions.options = [];
        amOptions.alpha   = [];
        amOptions.beta    = [];
    end

    %%% Default method
    if isempty(amOptions.method)
        amOptions.method = 'Standard';
    end

    if isempty(amOptions.options)
        amOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );    % default fsolve options
    end

    %%% Parameters definition
    switch amOptions.method
        case 'Standard'
            alpha = 24;             % standard alpha parameter
            beta = [9 19 -5 1];     % standard beta parameter

            if not(isempty(amOptions.alpha)) || not(isempty(amOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    amOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(amOptions.alpha) || isempty(amOptions.beta)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(amOptions.alpha),[1,1])) || not(isequal(size(amOptions.beta),[1,4]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = amOptions.alpha;    % custom alpha parameter
            beta = amOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-3)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), ...
                    zeros(dimSys,length(t)-3)];     % fval vector allocation
    feval = feval + 3*dimSys;                       % function evaluation counter update

    %%% AM4 loop
    for i = 3 : (length(t)-1)         % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        tp = t(i) + h;
        fp = @(xp) x(:,i) + h/alpha * (beta(1)*f(xp,tp) + beta(2)*fk1 + beta(3)*fk2 + beta(4)*fk3) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), amOptions.options);
        feval = feval + info.funcCount; 

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [amOptions] = amSettings(order,startup,method,options,alpha,beta)
%AM SETTINGS - Create the options struct for Adams Moulton methods
%
%   Syntax:
%       [amOptions] = amSettings(order,startup,method,options,alpha,beta)
%
%   Input:
%       order,          double:  select AB method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       options(*),     fsolve:  optimization options, for details see optimoptions.m
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,n]:  insert custom parameter beta for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AM method:
%           - amOptions.order = order
%           - amOptions.startup = startup;
%           - amOptions.method = method
%           - amOptions.options = options
%           - amOptions.alpha = alpha
%           - amOptions.beta = beta
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside am# functions
%       options: set to don't display iteration by default inside am# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%


    %%% Default value for optional input
    if nargin < 6    
        beta = [];
        if nargin < 5
            alpha = [];
            if nargin < 4
                options = [];
                if nargin < 3
                    method = [];
                    if nargin < 2
                        startup = [];
                    end
                end
            end
        end
    end

    %%% Struct definition
    amOptions = struct;
        amOptions.order   = order;
        amOptions.startup = startup;
        amOptions.method  = method;
        amOptions.options = options;
        amOptions.alpha   = alpha;
        amOptions.beta    = beta;

end

function [x,t,info] = backwardDifferenceFormula(f,x0,tmax,h,bdfOptions,visualConfig)
%BACKWARD DIFFERENCE FORMULA - Backward difference formula method selection and application
%
%   Syntax:
%       [x,t,info] = backwardDifferenceFormula(f,x0,tmax,h,bdfOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       bdfOtpions,     struct:  see bdfSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    %%% Dimension Check for initial guess
    if size(x0,2) > 1 && visualConfig == true
        % warning to remember that initial guess correct input format
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step\n',h);
    end

    if size(x0,2) > bdfOptions.order
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')

    elseif size(x0,2) < bdfOptions.order 
        if isempty(bdfOptions.startup)
            error('Need to define startup method to retrive a valid initial guess')
        end

    else
        infoStartup = [];   % initialization of the infoStartup variable used in info
        if not(isempty(bdfOptions.startup))
            warning('Startup method defined will not be used, x0 is already valid')
        end
    end

    %%% Adams Bashforth Moulton order selector based on 'bdfOptions'
    switch bdfOptions.order
        case 1
            [x,t,infoMethod] = bdf1(f,x0,tmax,h,bdfOptions);

        case 2
            % retrive a valid initial guess based on 'bdfOptions'
            if size(x0,2) < 2
                guessOrder = 2;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,bdfOptions);
            end
            
            [x,t,infoMethod] = bdf2(f,x0,tmax,h,bdfOptions);
            
        case 3
            % retrive a valid initial guess based on 'bdfOptions'
            if size(x0,2) < 3
                guessOrder = 3;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,bdfOptions);
            end

            [x,t,infoMethod] = bdf3(f,x0,tmax,h,bdfOptions);

        case 4
            % retrive a valid initial guess based on 'bdfOptions'
            if size(x0,2) < 4
                guessOrder = 4;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,bdfOptions);
            end

            [x,t,infoMethod] = bdf4(f,x0,tmax,h,bdfOptions);

        otherwise
            error('Please insert a valid order as input\n');
    end

    if visualConfig == true
        plot(t,x,'-');
    end

    if nargout == 3
        info = struct;
            info.fvalVec  = infoMethod.fvalVec; 
            info.implicit = true;

        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost  = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
        else
            info.timeCost  = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
        end
    end

end

function [x,t,info] = bdf1(f,x0,tmax,h,bdfOptions)
%BDF1 - Backward Difference Formula method of order 1
%
%   Syntax:
%       [x,t,info] = bdf1(f,x0,tmax,h,bdfOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       bdfOtpions(*),  struct:  see bdfSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       bdfOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        bdfOptions.method  = [];
        bdfOptions.options = [];
        bdfOptions.alpha   = [];
        bdfOptions.beta    = [];
        bdfOptions.gamma   = [];
    end

    %%% Default method
    if isempty(bdfOptions.method)
        bdfOptions.method = 'Standard';
    end

    if isempty(bdfOptions.options)
        bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );   % default fsolve options
    end

    %%% Parameters definition
    switch bdfOptions.method
        case 'Standard'
            alpha = 1;      % standard alpha parameter
            beta = 1;       % standard beta parameter
            gamma = 1;      % standard gamma parameter

            if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    bdfOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,1]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = bdfOptions.alpha;   % custom alpha parameter
            beta = bdfOptions.beta;     % custom beta parameter
            gamma = bdfOptions.gamma;   % custom gamma parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)]; % fval vector allocation
    feval = feval + dimSys;                     % function evaluation counter update

    %%% BDF1 loop
    for i = 1 : (length(t)-1)       % main loop of the method
        xk1 = x(:,i);
        tp = t(i) + h;
        fp = @(xp) 1/alpha * ( gamma(1)*xk1 + beta*h*f(xp,tp) ) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), bdfOptions.options);
        feval = feval + info.funcCount;

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [x,t,info] = bdf2(f,x0,tmax,h,bdfOptions)
%BDF2 - Backward Difference Formula method of order 2
%
%   Syntax:
%       [x,t,info] = bdf2(f,x0,tmax,h,bdfOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,2]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       bdfOtpions(*),  struct:  see bdfSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       bdfOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        bdfOptions.method  = [];
        bdfOptions.options = [];
        bdfOptions.alpha   = [];
        bdfOptions.beta    = [];
        bdfOptions.gamma   = [];
    end

    %%% Default method
    if isempty(bdfOptions.method)
        bdfOptions.method = 'Standard';
    end

    if isempty(bdfOptions.options)
        bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );   % default fsolve options
    end

    %%% Parameters definition
    switch bdfOptions.method
        case 'Standard'
            alpha = 3;          % standard alpha parameter
            beta = 2;           % standard beta parameter
            gamma = [4 -1];     % standard gamma parameter
            
            if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    bdfOptions.method);     % warning for parameters matrix unused
            end
            
        case 'Custom'
            if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)
                error('No parameters matrix has been given as input\n');  % missing parameters matrix
                
            elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,2]))
                error('Parameters matrix dimensions are invalid\n');      % parameters matrix with wrong size
                
            end
            
            alpha = bdfOptions.alpha;   % custom alpha parameter
            beta = bdfOptions.beta;     % custom beta parameter
            gamma = bdfOptions.gamma;   % custom gamma parameter
            
        otherwise
            error('Insert a valid method as input');
            
    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-2)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), ...
                    zeros(dimSys,length(t)-2)]; % fval vector allocation
    feval = feval + 2*dimSys;                   % function evaluation counter update

    %%% BDF2 loop
    for i = 2 : (length(t)-1)       % main loop of the method
        xk1 = x(:,i  );
        xk2 = x(:,i-1);
        tp = t(i) + h;
        fp = @(xp) 1/alpha * ( gamma(1)*xk1 + gamma(2)*xk2 + beta*h*f(xp,tp) ) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), bdfOptions.options);
        feval = feval + info.funcCount;

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [x,t,info] = bdf3(f,x0,tmax,h,bdfOptions)
%BDF3 - Backward Difference Formula method of order 3
%
%   Syntax:
%       [x,t,info] = bdf3(f,x0,tmax,h,bdfOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,3]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       bdfOtpions(*),  struct:  see bdfSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       bdfOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        bdfOptions.method  = [];
        bdfOptions.options = [];
        bdfOptions.alpha   = [];
        bdfOptions.beta    = [];
        bdfOptions.gamma   = [];
    end

    %%% Default method
    if isempty(bdfOptions.method)
        bdfOptions.method = 'Standard';
    end

    if isempty(bdfOptions.options)
        bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );   % default fsolve options
    end

    %%% Parameters definition
    switch bdfOptions.method
        case 'Standard'
            alpha = 11;             % standard alpha parameter
            beta = 6;               % standard beta parameter
            gamma = [18 -9 2];      % standard gamma parameter

            if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    bdfOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,3]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = bdfOptions.alpha;   % custom alpha parameter
            beta = bdfOptions.beta;     % custom beta parameter
            gamma = bdfOptions.gamma;   % custom gamma parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-3)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), ...
                    zeros(dimSys,length(t)-3)]; % fval vector allocation
    feval = feval + 3*dimSys;                   % function evaluation counter update

    %%% BDF3 loop
    for i = 3 : (length(t)-1)       % main loop of the method
        xk1 = x(:,i  );
        xk2 = x(:,i-1);
        xk3 = x(:,i-2);
        tp = t(i) + h;
        fp = @(xp) 1/alpha * ( gamma(1)*xk1 + gamma(2)*xk2 + gamma(3)*xk3 + beta*h*f(xp,tp) ) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), bdfOptions.options);
        feval = feval + info.funcCount;

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [x,t,info] = bdf4(f,x0,tmax,h,bdfOptions)
%BDF4 - Backward Difference Formula method of order 4
%
%   Syntax:
%       [x,t,info] = bdf4(f,x0,tmax,h,bdfOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,4]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       bdfOtpions(*),  struct:  see bdfSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       bdfOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        bdfOptions.method  = [];
        bdfOptions.options = [];
        bdfOptions.alpha   = [];
        bdfOptions.beta    = [];
        bdfOptions.gamma   = [];
    end

    %%% Default method
    if isempty(bdfOptions.method)
        bdfOptions.method = 'Standard';
    end

    if isempty(bdfOptions.options)
        bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );   % default fsolve options
    end

    %%% Parameters definition
    switch bdfOptions.method
        case 'Standard'
            alpha = 25;             % standard alpha parameter
            beta = 12;              % standard beta parameter
            gamma = [48 -36 16 3];  % standard gamma parameter

            if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    bdfOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,4]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = bdfOptions.alpha;   % custom alpha parameter
            beta = bdfOptions.beta;     % custom beta parameter
            gamma = bdfOptions.gamma;   % custom gamma parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-4)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), f(x(:,4),t(4)), ...
                    zeros(dimSys,length(t)-4)]; % fval vector allocation
    feval = feval + 4*dimSys;                   % function evaluation counter update

    %%% AB4 loop
    for i = 4 : (length(t)-1)       % main loop of the method
        xk1 = x(:,i  );
        xk2 = x(:,i-1);
        xk3 = x(:,i-2);
        xk4 = x(:,i-3);
        tp = t(i) + h;
        fp = @(xp) 1/alpha * ( gamma(1)*xk1 + gamma(2)*xk2 + gamma(3)*xk3 + gamma(4)*xk4 + beta*h*f(xp,tp) ) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), bdfOptions.options);
        feval = feval + info.funcCount;

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

end

function [bdfOptions] = bdfSettings(order,startup,method,options,alpha,beta,gamma)
%BDF SETTINGS - Create the options struct for Backward Difference methods
%
%   Syntax:
%       [bdfOptions] = bdfSettings(order,startup,method,options,alpha,beta,gamma)
%
%   Input:
%       order,          double:  select AB method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       options(*),     fsolve:  optimization options, for details see optimoptions.m
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,1]:  insert custom parameter beta for AB method
%       gamma(*),  double[1,n]:  insert custom parameter gamma for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AB method:
%           - amOptions.order = order
%           - amOptions.startup = startup;
%           - amOptions.method = method
%           - amOptions.options = options
%           - amOptions.alpha = alpha
%           - amOptions.beta = beta
%           - amOptions.gamma = gamma
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside ab# functions
%       options: set to don't display iteration by default inside am# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%       gamma:   set as empty by default
%


    %%% Default value for optional input
    if nargin < 7
        gamma = [];
        if nargin < 6    
            beta = [];
            if nargin < 5
                alpha = [];
                if nargin < 4
                    options = [];
                    if nargin < 3
                        method = [];
                        if nargin < 2
                            startup = [];
                        end
                    end
                end
            end
        end
    end

    %%% Struct definition
    bdfOptions = struct;
        bdfOptions.order   = order;
        bdfOptions.startup = startup;
        bdfOptions.method  = method;
        bdfOptions.options = options;
        bdfOptions.alpha   = alpha;
        bdfOptions.beta    = beta;
        bdfOptions.gamma   = gamma;

end

function [cmap] = graphicSettings()
%GRAPHIC SETTINGS - set grapghics options for the visual output
%   

    % Graphic settings
    set(0,'defaulttextinterpreter','latex');  
    set(0,'defaultAxesTickLabelInterpreter','latex');  
    set(0,'defaultLegendInterpreter','latex');
    set(0,'DefaultLineLineWidth',1.2);
    warning('off','MATLAB:fplot:NotVectorized');

    % Color map definition
    cmap = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
    
    % Create folder for figure export
    if ~exist('figure','dir')
        mkdir('figure');
    end
    
end

function [x,t,info] = iex4(f,x0,tmax,h,iOptions,visualConfig)
%IEX4 - Implicit Extrapolation method application
%
%   Syntax:
%       [x,t,info] = iex4(f,x0,tmax,h,iOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       iOtpions(*),    struct:  see iSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       iOtpions:      set with default solver 'fsolve'
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    if nargin < 5
        iOptions.method  = [];
        iOptions.options = [];
    end

    %%% Default submethod
    if isempty(iOptions.method)
        iOptions.method = 'fsolve';      % default implicit solver
    end

    if isempty(iOptions.options)
        iOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );  % default fsolve options
    end

    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected')
    end
                    
    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    alpha = [-1/6, 4, -27/2, 32/3];                 % iex4 parameters definition
    x = [x0 zeros(dimSys,length(t)-1)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];     % fval vector allocation
    feval = feval + dimSys;                         % function evaluation counter update

    %%% IEX4 loop
    for i = 1 : (length(t)-1)         % main loop of the method

        %1st predictor
        tk1 = t(i) + h;
        fk1 = @(k1) x(:,i) + h * f(k1,tk1) - k1;
        [xp1,~,conv(1),info] = fsolve(fk1, x(:,i), iOptions.options);         % initial guess xk
        feval = feval + info.funcCount;

        %2nd predictor
        tk2a = t(i) + h/2;
        fk2a = @(k2a) x(:,i) + h/2 * f(k2a,tk2a) - k2a;
        [k2a,~,conv(2),info] = fsolve(fk2a, x(:,i), iOptions.options);        % initial guess xk
        feval = feval + info.funcCount;

        tk2 = tk1;
        fk2 = @(k2) k2a + h/2 * f(k2,tk2) - k2;
        [xp2,~,conv(3),info] = fsolve(fk2, k2a, iOptions.options);            % initial guess k2a
        feval = feval + info.funcCount;

        %3rd predictor
        tk3a = t(i) + h/3;
        fk3a = @(k3a) x(:,i) + h/3 * f(k3a,tk3a) - k3a;
        [k3a,~,conv(4),info] = fsolve(fk3a, x(:,i), iOptions.options);        % initial guess xk
        feval = feval + info.funcCount;

        tk3b = t(i) + h * 2/3;
        fk3b = @(k3b) k3a + h/3 * f(k3b,tk3b) - k3b;
        [k3b,~,conv(5),info] = fsolve(fk3b, k3a, iOptions.options);           % initial guess k3a
        feval = feval + info.funcCount;

        tk3 = tk1;
        fk3 = @(k3) k3b + h/3 * f(k3,tk3) - k3;
        [xp3,~,conv(6),info] = fsolve(fk3, k3b, iOptions.options);            % initial guess k3b
        feval = feval + info.funcCount;

        %4th predictor
        tk4a = t(i) + h/4;
        fk4a = @(k4a) x(:,i) + h/4 * f(k4a,tk4a) - k4a;
        [k4a,~,conv(7),info] = fsolve(fk4a, x(:,i), iOptions.options);        % initial guess xk
        feval = feval + info.funcCount;

        tk4b = t(i) + h * 2/4;
        fk4b = @(k4b) k4a + h/4 * f(k4b,tk4b) - k4b;
        [k4b,~,conv(8),info] = fsolve(fk4b, k4a, iOptions.options);           % initial guess k4a
        feval = feval + info.funcCount;

        tk4c = t(i) + h * 3/4;
        fk4c = @(k4c) k4b + h/4 * f(k4c,tk4c) - k4c;
        [k4c,~,conv(9),info] = fsolve(fk4c, k4b, iOptions.options);           % initial guess k4b
        feval = feval + info.funcCount;

        tk4 = tk1;
        fk4 = @(k4) k4c + h/4 * f(k4,tk4) - k4;
        [xp4,~,conv(10),info] = fsolve(fk4, k4c, iOptions.options);            % initial guess k4c
        feval = feval + info.funcCount;

        % corrector
        x(:,i+1) = alpha(1) * xp1 + alpha(2) * xp2 + alpha(3) * xp3 + alpha(4) * xp4;
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;     % function evaluation counter update

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end

    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end

    %%% plot
    if visualConfig == true
        plot(t,x,'o-');
    end

end

function [iOptions] = iSettings(method,options)
%I SETTINGS - Create the options struct for iex4 methods
%
%   Syntax:
%       [iOptions] = iSettings(method,options)
%
%   Input:
%       method(*),        char:  define the solver (available 'fsolve' only in the current realese)
%       options(*),     fsolve:  optimization options, for details see optimoptions.m
%
%   Output:
%       iOptions,  struct:  contains settings for IEX4 method:
%           - iOptions.method = method
%           - iOptions.options = options
%
%   Default settings for optional input (*):
%       method:  set as 'fsolve' by default inside iex4 function
%       options: set to don't display iteration by default inside iex4 function
%


    %%% Default value for optional input
    if nargin < 2
        options = [];
        if nargin < 1
            method = [];
        end
    end

    %%% Struct definition
    iOptions = struct;
        iOptions.method  = method;
        iOptions.options = options;

end

function [solution,converge,info] = newton(f,xGuess,method,toll,nmax,iter,visualConfig)
%NEWTON - Compute zeros of given function using the Newton's method
%
%   Syntax:
%       [solution,converge,info] = newton(f,xGuess,method,toll,nmax,visualConfig)
%
%   Input:
%       f,   function(x1,...,xn):  function for which to compute zeros
%       xGuess,      double[n,1]:  generic initial guess 
%       method(*),          char:  chose between:
%               - 's' compute jacobian using symbolic math
%               - 'f' compute jacobian using forward finite difference
%               - 'c' compute jacobian using centered finite difference
%       toll(*),          double:  tollerance on the error between iterations
%       nmax(*),          double:  maximum number of iterations
%       iter(*),          double:  number of iteration to compute time cost
%       visualConfig(*),  struct:  set visual output of the function
%               - visualConfig.print: display solutions and some info
%               - visualConfig.plot: plot the convergence plot
%
%   Output:
%       solution,           double:  solution vector
%       convergence,   double[1,m]:  time istant associated to solutions
%       info,               struct:  information on method used:
%           - info.iteration,           double:  number of iteration (m)
%           - info.solutionVector, double[n,m]:  solution vector for each iteration
%           - info.errorVector,    double[1,m]:  relative error between iterations
%           - info.timeCost,            double:  time spent
%           - info.fevalCost,           double:  # of function evaluations
%           - info.methodUsed,            char:  string describing the method used
%           - info.absError,            double:  absolute error on mathematical zero
%
%   Default settings for optional input (*):
%       method:        set as 'f' by default
%       toll:          set as 1e-12 by default
%       nmax:          set as 1e3 by default
%       visualConfig:  set as true by default
%


    %%% Prelimirary check
    if nargin < 6                       % default options without visual config flags
        visualConfig.print = true;
        visualConfig.plot = true;
    end

    if nargin < 5 || isempty(nmax)      % default max number of iteration
        nmax = 1e3;
    end

    if nargin < 4 || isempty(toll)      % default tollerance
        toll = 1e-12;
    end

    if nargin < 5 || isempty(iter)      % default tollerance
        iter = 10;
    end

    if method ~= 's' && method ~= 'f' && method ~= 'c'
        error ('Method not valid, please insert a valid method as input')
    end

    timerStart = tic;               % timer start
    for timerIter = 1:iter

        %%% Initialization
        feval = 0;                      % function evaluation counter starts
        dimSys = length(xGuess);        % function evaluation step
        err = toll + 1;                 % error
        k = 0;                          % iterator
        errVect = [];                   % relative error vector
        xkVect = [];                    % solutions vector
        xGuessCell = num2cell(xGuess);  % inital guess cell array

        %%% Newton method body
        % Check if initial guess is already a solution   
        if f(xGuessCell{:,1}) == zeros(length(xGuess),1)
            elapsedTime = toc(timerStart);
            solution = xGuess;
            converge = true;

            if nargout == 3
                info = struct;
                info.iteration      = 0;
                info.solutionVector = solution;
                info.errorVector    = f(xGuessCell{:,1});
                info.timeCost       = elapsedTime;
                info.fevalCost      = dimSys;
                info.methodUsed     = 'none';
            end

            % output print on command window
            if visualConfig.print == true
                fprintf('Function zeros are exactly equal to the initial guess\n');
                solPrint = sprintf(' %.4f\n',solution);
                fprintf('\n Zeros computed:\n%s \n',solPrint)
            end

            % function in this "lucky" case ends here
            return

        else

            feval = feval + dimSys;     % function evaluation used to check if the guess is already a solution
            xk = xGuess;                % start point of the iteration process

            % symbolic manipulation to compute jacobian matrix
            if method == 's'
                fSym = sym(f);
                JSym = jacobian(fSym);
                J = matlabFunction(JSym);
            end

            % Iterative loop to find the solution starts
            while k < nmax && err > toll

                switch method
                    case 's'
                        % Symbolic math method
                        xkCell = num2cell(xk);          % xk convertion in cell array
                        fk = f(xkCell{:,1});            % function f evaluated in xk
                        feval = feval + dimSys;         % function evaluation counter update

                        Jk = J(xkCell{:,1});                % jacobian f evaluated in xk
                        feval = feval + dimSys^(dimSys);    % function evaluation counter update

                        methodString = 'SYMBOLIC MATH';  % method used saved for the output

                    case 'f'
                        % Finite forward difference method
                        epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));     % epsilon coputed at step xk (vector: dimSys,1)
                        Jk = zeros(dimSys);                             % jacobian initialization

                        xkCell = num2cell(xk);                          % xk convertion in cell array
                        fk = f(xkCell{:,1});                            % function f evaluated in xk
                        feval = feval + dimSys;                         % function evaluation counter update

                        % loop to evaluate jacobian's i^th column in xk using forward difference
                        for i=1:dimSys
                            epsilonVect = (zeros(dimSys,1));            % epsilon vector initialization
                            epsilonVect(i,1) = epsilon(i,1);            % epsilon vector definition

                            xkForward =  xk + epsilonVect;              % xk+epsilon computation in vector
                            xkForwardCell = num2cell(xkForward);        % xk+epsilon convertion in cell array
                            fkForward = f(xkForwardCell{:,1});          % function f evaluated in xk+epsilon
                            feval = feval + dimSys;                     % function evaluation counter update

                            Jk(:,i) = (fkForward - fk) / epsilon(i,1);  % i^th columns of jacobian evaluated in xk
                        end

                        methodString = 'FORWARD DIFF';                  % method used saved for the output

                    case 'c'
                        % Finite centered difference method
                        epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));     % epsilon coputed at step xk (vector: dimSys,1)
                        Jk = zeros(dimSys);                             % jacobian initialization

                        xkCell = num2cell(xk);                          % xk convertion in cell array
                        fk = f(xkCell{:,1});                            % function f evaluated in xk
                        feval = feval + dimSys;                         % function evaluation counter update

                        % loop to evaluate jacobian's i^th column in xk using centered difference
                        for i=1:dimSys
                            epsilonVect = (zeros(dimSys,1));            % epsilon vector initialization
                            epsilonVect(i,1) = epsilon(i,1);            % epsilon vector definition

                            xkForward =  xk + epsilonVect;              % xk+epsilon computation in vector
                            xkForwardCell = num2cell(xkForward);        % xk+epsilon convertion in cell array
                            fkForward = f(xkForwardCell{:,1});          % function f evaluated in xk+epsilon
                            feval = feval + dimSys;                     % function evaluation counter update

                            xkBackward =  xk - epsilonVect;             % xk-epsilon computation in vector
                            xkBackwardCell = num2cell(xkBackward);      % xk-epsilon convertion in cell array
                            fkBackward = f(xkBackwardCell{:,1});        % function f evaluated in xk-epsilon
                            feval = feval + dimSys;                     % function evaluation counter update

                            Jk(:,i) = (fkForward - fkBackward) / (2*epsilon(i,1));   % i^th columns of jacobian evaluated in xk
                        end

                        methodString = 'CENTERED DIFF';  % method used saved for the output

                end

                % check if f is sufficiently well-behaved
                if rank(Jk) == zeros(length(xGuess))
                    error('Jk has become zero, f not sufficiently well-behaved near solution\n');

                else
                    y = Jk \ (-fk);
                    xNew = xk + y;              % solutions for the k^th iteration
                    err = norm(y);              % error computation at k^th iteration

                    errVect = [errVect, err];   % relative error vector update
                    xkVect  = [xkVect, xNew];   % solution vector update
                    k = k + 1;                  % iteration index update
                    xk = xNew;                  % starting point update for the next iteration

                end
            end
        end
        absError = norm(f(xk(1),xk(2)));        % absolute final error
    end
    elapsedTime = toc(timerStart);          % timer stop

    %%% Output organization
    solution = xkVect(:,end);   % final result
    if err > toll           % convergence flag
        converge = false;
    else
        converge = true;
    end

    if nargout == 3
        info = struct;      % info struct for detailed output
            info.iteration      = k;
            info.solutionVector = xkVect;
            info.errorVector    = errVect;
            info.timeCost       = elapsedTime/iter;
            info.fevalCost      = feval;
            info.methodUsed     = methodString;
            info.absError       = absError;
    end

    % output print on command window
    if visualConfig.print == true
        if converge == true
            fprintf('Function zeros has been successfully computed using %s method \n',methodString);
            fprintf('After %d interations the error %.3d is less than the tollerance %.3d \n', ...
                k, absError, toll);
        else
            fprintf('Function zeros has NOT been computed using %s method \n',methodString);
            fprintf('After %d interations the error %.3d is still grater than the tollerance %.3d \n', ...
                k, absError, toll);
        end
        solPrint = sprintf(' %.4f\n',solution);
        fprintf('\n Zeros computed:\n%s \n',solPrint)
    end

    % convergence plot
    if visualConfig.plot == true
        iterationVector = (1:k);
        semilogy(iterationVector,errVect,'o-')
    end

end

function [x,t,info] = rk1(f,x0,tmax,h,rkOptions)
%RK1 - Runge-Kutta method of order 1
%
%   Syntax:
%       [x,t,info] = rk1(f,x0,tmax,h,rkOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions(*),   struct:  see rkSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       rkOptions:  set with default 'FowardEuler' alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end

    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'FowardEuler';   % FowardEuler parameters sets as default
    end

    if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
        warning('Parameters matrix unused, standard %s parameters are used instead',...
            rkOptions.method);      % warning for parameters matrix unused
    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];           % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% RK1 loop
    for i=1:(length(t)-1)
        x(:,i+1) = x(:,i) + h * f(x(:,i),t(i));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;      % function evaluation counter update
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;               % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = rk2(f,x0,tmax,h,rkOptions)
%RK2 - Runge-Kutta method of order 2
%
%   Syntax:
%       [x,t,info] = rk2(f,x0,tmax,h,rkOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions(*),   struct:  see rkSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       rkOptions:  set with default 'Heun' alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end

    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'Heun';          % Heun parameters sets as default
    end

    %%% Parameters definition
    switch rkOptions.method
        case 'Heun'
            alpha2 = [1 1]';                 % rk alpha matrix
            beta2 = [1, 0; 0.5, 0.5];        % rk beta matrix

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'MidPoint'
            alpha2 = [0.5 1]';               % rk alpha matrix
            beta2 = [0.5, 0; 0, 1];          % rk beta matrix

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(rkOptions.alpha) || isempty(rkOptions.beta)
                error('No parameters matrix has been given as input')   % missing parameters matrix

            elseif not(isequal(size(rkOptions.alpha),[2,1])) || not(isequal(size(rkOptions.beta),[2,2]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size
                
            end

            alpha2 = rkOptions.alpha;        % custom alpha parameter
            beta2 = rkOptions.beta;          % custom beta parameter

            %%% Check if the custom parameters are valid
            eC1 = beta2(2,1) + beta2(2,2) == 1;
            eC2 = 2 * beta2(1,1) * beta2(2,2) == 1;
            eC3 = 2 * alpha2(1,1) * beta2(2,2) == 1;

            if eC1 ~= true ||  eC2 ~= true ||  eC3 ~= true
                error('Insert a valid parameters matrix as input');
            end

        otherwise
            error('Insert a valid method as input');
    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x =  [x0 zeros(dimSys,length(t)-1)];          % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% RK2 loop
    for i=1:(length(t)-1)
        fk = fvalVec(:,i);
        xp = x(:,i) + beta2(1,1) * h * fk;    % x(i)=xk && t(i)=tk
        tp = t(i) + alpha2(1,1) * h;
        x(:,i+1) = x(:,i) + alpha2(2,1) * h * (beta2(2,1) * fk + beta2(2,2) * f(xp,tp));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;             % function evaluation counter update
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = rk3(f,x0,tmax,h,rkOptions)
%RK3 - Runge-Kutta method of order 3
%
%   Syntax:
%       [x,t,info] = rk3(f,x0,tmax,h,rkOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions(*),   struct:  see rkSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       rkOptions:  set with default alpha and beta
%

    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end

    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'Default';           % Default parameters sets as default
    end

    %%% Parameters definition
    switch rkOptions.method
        case 'Default'
            alpha3 = [1/3 2/3 1]';           
            beta3 = diag([1/3 2/3 3/4 ]);       
            beta3(3,:) = [1/4 0 3/4];

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(rkOptions.alpha) || isempty(rkOptions.beta)            
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(rkOptions.alpha),[3,1])) || not(isequal(size(rkOptions.beta),[3,3]))            
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha3 = rkOptions.alpha;        % custom alpha parameter
            beta3 = rkOptions.beta;          % custom beta parameter

        otherwise
            error('Insert a valid method as input');
    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x =  [x0 zeros(dimSys,length(t)-1)];          % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% RK3 loop
    for i=1:(length(t)-1)                              % calculation loop
        fk = fvalVec(:,i);
        xp1 = x(:,i) + beta3(1,1) * h * fk;            % x(i)=xk && t(i)=tk
        tp1 = t(i) + alpha3(1,1) * h;

        xp2 = x(:,i) + beta3(2,2) * h * f(xp1,tp1);    % x(i)=xk && t(i)=tk
        tp2 = t(i) + alpha3(2,1) * h;

        x(:,i+1) = x(:,i) + alpha3(3,1) * h * ( beta3(3,1) * fk ...
            + beta3(3,2) * f(xp1,tp1) ...
            + beta3(3,3) * f(xp2,tp2));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 3*dimSys;             % function evaluation counter update
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [x,t,info] = rk4(f,x0,tmax,h,rkOptions)
%RK4 - Runge-Kutta method of order 4
%
%   Syntax:
%       [x,t,info] = rk4(f,x0,tmax,h,rkOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions(*),   struct:  see rkSettings.m for details 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       rkOptions:  set with default 'Runge-Kutta' alpha and beta
%

    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end

    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'Runge-Kutta';      % Runge-Kutta parameters sets as default
    end

    %%% Parameters definition
    switch rkOptions.method
        case 'Runge-Kutta'
            alpha4 = [0.5 0.5 1 1]';           % rk alpha matrix
            beta4 = diag([0.5 0.5 1 1/6]);     % rk beta matrix
            beta4(4,:) = [1/6 1/3 1/3 1/6];

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(rkOptions.alpha) || isempty(rkOptions.beta)            
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(rkOptions.alpha),[4,1])) || not(isequal(size(rkOptions.beta),[4,4]))            
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha4 = rkOptions.alpha;        % custom alpha parameter
            beta4 = rkOptions.beta;          % custom beta parameter

        otherwise
            error('Insert a valid method as input');
    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x =  [x0 zeros(dimSys,length(t)-1)];          % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% RK4 loop
    for i=1:(length(t)-1)                              % calculation loop
        fk = fvalVec(:,i);
        xp1 = x(:,i) + beta4(1,1) * h * fk;            % x(i)=xk && t(i)=tk
        tp1 = t(i) + alpha4(1,1) * h;

        xp2 = x(:,i) + beta4(2,2) * h * f(xp1,tp1);    % x(i)=xk && t(i)=tk
        tp2 = t(i) + alpha4(2,1) * h;

        xp3 = x(:,i) + beta4(3,3) * h * f(xp2,tp2);    % x(i)=xk && t(i)=tk
        tp3 = t(i) + alpha4(3,1) * h;

        x(:,i+1) = x(:,i) + alpha4(4,1) * h * ( beta4(4,1) * fk ...
            + beta4(4,2) * f(xp1,tp1) ...
            + beta4(4,3) * f(xp2,tp2)...
            + beta4(4,4) * f(xp3,tp3) );
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 4*dimSys;             % function evaluation counter update
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end

function [rkOptions] = rkSettings(order,method,alpha,beta,iterations)
%RK SETTINGS - Create the options struct for Runge-Kutta methods
%
%   Syntax:
%       [rkOptions] = rkSettings(order,method,alpha,beta)
%
%   Input:
%       order,          double:  select RK method order
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       alpha(*),  double[1,1]:  insert custom parameter alpha for RK method
%       beta(*),   double[1,n]:  insert custom parameter beta for RK method
%       iterations(*),  double:  number of iterations to compute average CPU time
%
%   Output:
%       rkOptions,  struct:  contains settings for RK method:
%           - rkOptions.order = order
%           - rkOptions.method = method
%           - rkOptions.alpha = alpha
%           - rkOptions.beta = beta
%           - rkOptions.iterations = iterations
%
%   Default settings for optional input (*):
%       method:  set as 'Standard' by default inside rk# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%       iterations: set as 1 by default
%

    %%% Default value for optional input
    if nargin < 5
        iterations = 1;
        if nargin < 4
            beta = [];
            if nargin < 3
                alpha = [];
                if nargin < 2
                    method = [];
                end
            end
        end
    end

    %%% Struct definition
    rkOptions = struct;
        rkOptions.order  = order;
        rkOptions.method = method;
        rkOptions.alpha  = alpha;
        rkOptions.beta   = beta;
        rkOptions.iterations = iterations;

end

function [x,t,info] = rungeKutta(f,x0,tmax,h,rkOptions,visualConfig)
%RUNGE KUTTA - Runge-Kutta method selection and application
%
%   Syntax:
%       [x,t,info] = rungeKutta(f,x0,tmax,h,rkOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions,      struct:  see rkSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%           - info.iterations,    double:
%           - info.avgTimeCost,  double:  averge time spent on 'avgTime' iterations
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
    end
    
    %%% Runge-Kutta method selector based on 'rkOptions'
    switch rkOptions.order
        case 1 
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk1(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        case 2
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk2(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        case 3
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk3(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        case 4
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk4(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        otherwise
            error('Please insert a valid method as input\n');
    end
    info.iterations = rkOptions.iterations;

    if visualConfig == true
        plot(t,x,'-');
    end

end

function [F,guess,degVec,hVec] = selectOp(method,dimSys,params)
%SELECT OP - Finds the F(h,alpha) operator and associated settings usefull
% to compute stability region of the desired integration method 
%
%    Syntax:
%       [F,guess,degVec,hVec] = selectOp(method,dimSys,params)
%
%   Input:
%       method,       char:  select method to obtain the operator
%       dimSys(*),  double:  analyzed system dimension 
%       params(*),  double:  input used in BI2 method (sets theta value)
%
%   Output:
%       F,      function(h,A):  corresponding method operator
%       guess,         double:  eductaed initial guess used for fzero in stabRegion.m
%       degVec,   double[1,2]:  starting and ending alpha used in stability region analysis
%       hVec,     double[1,2]:  starting and ending h used in stability guess analysis
%
%   Default settings for optional input (*):
%       dimSys:     set to 2 by default for stability region analysis
%       params:     set to 0.4 by default
%
%   Note:   'guess','degVec','hVec' has been pre defined throught an iterarive
%           based on graphic analysis obtained using the function stabGuess.m
%


    %%% Optional input definition
    if nargin < 2
        dimSys = 2;
        if nargin < 3       
            params = 0.4;
        end
    end
                 
    %%% Operator selection
    I = eye(dimSys);
    switch method
        case 'RK1'
            F = @(h,A) I + (h*A);
            guess = 2;
            degVec = [180 0];
            hVec = [0 5];
    
        case 'RK2'
            F = @(h,A) I + (h*A) + 1/2*(h*A)^2;
            guess = 2;
            degVec = [180 0];
            hVec = [0 5];
    
        case 'RK3'
            F = @(h,A) I + (h*A) + 1/2*(h*A)^2 + 1/6*(h*A)^3;
            guess = 2;
            degVec = [180 0];
            hVec = [0 5];
    
        case 'RK4'
            F = @(h,A) I + (h*A) + 0.5*(h*A)^2 + 1/6*(h*A)^3 + 1/24*(h*A)^4;
            guess = 2;
            degVec = [180 0];
            hVec = [0 5];
    
        case 'BI2'
            F = @(h,A) (I-(1-params)*h*A+((1-params)*h*A)^2/2) \ (I+params*h*A+(params*h*A)^2/2);
    
            if params < 0.5
                guess  = 10;
                degVec = [0 180];
                hVec = [0 12];
    
            elseif params >= 0.5
                guess  = 10;
                degVec = [180 0];
                hVec = [0 12];
            end
            
        case 'IEX4'
            F  = @(h,A) -1/6*((I-h*A/1)\I)^1 + 4*((I-h*A/2)\I)^2 - 27/2*((I-h*A/3)\I)^3 + 32/3*((I-h*A/4)\I)^4;
            guess  = 10;
            degVec = [0 180];
            hVec = [10 15];
    
        otherwise
            error('Insert a valid method')
    
    end

end

function [] = stabGuess(F,degStart,hVec)
%STAB GUESS - Plot the Stability problem function in the starting alpha
%
%   Syntax:
%       [] = stabGuess(F,degStart,hVec)
%
%   Input:
%       F,   function(h,A):  corresponding method operator
%       degStart,   double:  starting point of the stability analysis
%       hVec,  double[1,2]:  starting and ending h used in fplot
%
%     


    radStart = deg2rad(degStart);
    A = [0, 1; -1, 2*cos(radStart)];
    S = @(h) max(abs(eig(F(h,A)))) - 1;     % Stability problem definition
    
    fplot(S,hVec, Color = 'k',LineWidth=1.2);               % Initial guess check for stability region
    hold on;    grid on;    axis padded;    box on; 
    ylim([-1 1]);     yline(0, LineStyle = '--')
    xlabel('$h$');    ylabel('$S(h)$')                      % where S is the stability function
    tstring = sprintf('Stability function in %.0f',degStart);
    title (tstring);
    drawnow

end

function [xF,yF,hvec,stabEdge] = stabRegion(F,stepNum,degVec,guess,visualConfig)
%STAB REGION - Finds the stability region associated to the operator F(h,A)
%
%    Syntax:
%       [xF,yF,hvec,stabEdge] = stabRegion(F,stepNum,degVec,guess,visualConfig)
%
%   Input:
%       F,       function(h,A):  method operator to compute stability region
%       stepNum,        double:  number of intermidiate step between alpha limit 
%       degVec,    double[1,2]:  starting and ending alpha to compute stability region
%       guess,          double:  eductaed initial guess to find stability region
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       xF,       double[1, n]:  x coordinates of the h*lambda plane points between degVec values
%       yF,       double[1, n]:  y coordinates of the h*lambda plane points between degVec values
%       hVec,     double[1, n]:  h computed for the stability problem between degVec values
%       stabEdge, double[2,2n]:  stability region contour coordinates in all the h*lambda plane
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%

    
    %%% Optional input definition
    if nargin < 5
        visualConfig = true;
    end
    
    %%% Variables initialization
    Afun = @(a) [0, 1; -1, 2*cos(a)];
    alphaLim = deg2rad(degVec);
    alphaVec = linspace(alphaLim(1),alphaLim(2),stepNum);
    hvec = zeros(1,length(alphaVec));
    xF   = zeros(1,length(alphaVec));
    yF   = zeros(1,length(alphaVec));
    
    %%% Stability analysis loop
    for i=1:length(alphaVec)        
        
        alpha = alphaVec(i);
        A = Afun(alpha);                        % Local A matrix value
        S = @(h) max(abs(eig(F(h,A)))) - 1;     % Stability problem definition
        [hvec(i),~,conv]=fzero(S, guess);       % Stability problem resolution
              
        guess = hvec(i);        % Educated gusess based on the current solution       
        eig_iterVec = eig(A);   % Eigenvalue of continuos system
        xF(:,i) = hvec(i)*real(eig_iterVec(1)); % h*lambda x coordinates
        yF(:,i) = hvec(i)*imag(eig_iterVec(1)); % h*lambda y coordinates
        
        % convergence check
        if conv < 1
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end
    
    % Entire stability region coordinate vectors
    stabEdge = [xF,  flip(xF,2);yF, -flip(yF,2)];
    
    if visualConfig == true
        axis equal;     grid on;    box on;     hold on
        ax = gca;       ax.Layer = 'top';
            
        % final stability region computed plot
        plot(stabEdge(1,:), stabEdge(2,:),'k')
        if degVec(1) ~= 0       % explicit method case
            fill(stabEdge(1,:), stabEdge(2,:),[0.80 0.80 0.80])
        elseif degVec(1) == 0   % implicit method case
            ax.Color = [0.80 0.80 0.80];
            fill([-inf -inf inf inf],[-inf inf inf -inf],[0.80 0.80 0.80])
            fill(stabEdge(1,:), stabEdge(2,:),'w')
            grid on;
        end
        xline(0,'--');                  yline(0,'--')
        xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
        title ('Stability region')
    
    end

end

function [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,methodOptions)
%STARTUP GUESS - Solve the startup problem for multistep method
%
%   Syntax:
%       [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,methodOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       guessOrder,     double:  number of element in initial guess
%       h,              double:  time step of the integration
%       methodOptions,  struct:  see corresponding methodSettings.m for details
%
%   Output:
%       x,     double[n,m]:  solution vector
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%


    %%% Initialization
    if nargout == 2
        infoStartup = struct;
            infoStartup.timeCost  = [];
            infoStartup.fevalCost = [];
            infoStartup.implicit  = [];
    end
    visualConf = false;     % disable plot in the startup problem

    %%% Startup solver selector vased on 'methodOptions'
    switch methodOptions.startup
        case 'RK'
            % select the runge kutta associated to same order as multistep method
            rkStarter = rkSettings(guessOrder);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);
            
        case 'RK1'
            rkStarter = rkSettings(1);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

        case 'RK2'
            rkStarter = rkSettings(2);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

        case 'RK3'
            rkStarter = rkSettings(3);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);
            
        case 'RK4'
            rkStarter = rkSettings(4);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

        case 'IEX4'
            iStarter = iSettings();
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;      
            [x0,~,infoStartup] = iex4(f,x0,tmaxSingleStep,h,iStarter,visualConf);

        case 'AB'
            while size(x0,2) ~= guessOrder
                abStarter = abSettings(size(x0,2));
                tmaxLoop = abStarter.order * h;
                [x0,~,infoStartupLoop] = adamsBashforth(f,x0,tmaxLoop,h,abStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end

        case 'AM'
            while size(x0,2) ~= guessOrder
                amStarter = amSettings(size(x0,2));
                tmaxLoop = amStarter.order * h;
                [x0,~,infoStartupLoop] = adamsMoulton(f,x0,tmaxLoop,h,amStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end
        
        case 'ABM'
            while size(x0,2) ~= guessOrder
                abmStarter = abmSettings(size(x0,2));
                tmaxLoop = abmStarter.order * h;
                [x0,~,infoStartupLoop] = adamsBashforthMoulton(f,x0,tmaxLoop,h,abmStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end
        
        case 'BDF'
            while size(x0,2) ~= guessOrder
                bdfStarter = bdfSettings(size(x0,2));
                tmaxLoop = bdfStarter.order * h;
                [x0,~,infoStartupLoop] = backwardDifferenceFormula(f,x0,tmaxLoop,h,bdfStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end

        otherwise
            error('Please insert a valid method as input');
    end

end

function [Z1mesh,Z2mesh,t] = zerosGuess(f1,f2,meshMin,meshMax,meshStep,visualConfig)
%ZEROS GUESS - plot the two functions on a defined mesh to find an educated
% guess for 2x2 Newton method
%
%   Syntax:
%       [Z1mesh,Z2mesh,t] = zerosGuess(f1,f2,meshMin,meshMax,visualConfig)
%
%   Input:
%       f1,  function(x1,...,xn):  function 1 for which to compute zeros
%       f2,  function(x1,...,xn):  function 2 for which to compute zeros
%       meshMin(*),       double:  mesh lower limit
%       meshMax(*),       double:  mesh upper limit
%       meshStep(*),      double:  mesh step between to values
%       visualConfig(*),    bool:  set visual output of the function
%
%   Output:
%       Z1mesh,   double[n,m]:  mesh evaluated on f1
%       Z2mesh,   double[n,m]:  mesh evaluated on f2
%       t,   TiledChartLayout:  visual output object
%
%   Default settings for optional input (*):
%       meshMin:       set equal to -5 by default
%       meshMax:       set equal to -5 by default
%       meshStep:      set equal to 0.1 by defaylt
%       visualConfig:  set as true by default
%


    %%% Default value for optional input
    if nargin < 6
        visualConfig = true;
        if nargin < 5
            meshStep = 0.1;
            if nargin < 4
                meshMax = 5;
                if nargin < 3
                    meshMin = -5;
                end
            end
        end
    end
    
    %%% mesh definition and evaluation
    [X1mesh,X2mesh] = meshgrid(meshMin:meshStep:meshMax);
    Z1mesh = f1(X1mesh,X2mesh);
    Z2mesh = f2(X1mesh,X2mesh);    
    
    if visualConfig == true
        fig = figure();
            fig.Name = 'Function plot';
            %fig.Position = [1 1 4000 2000];    % to find values
        t = tiledlayout(1,2);
        tString = title(t,'Initial guess graphic search');
        tString.Interpreter = 'latex';
        
        %%% plot 3d of the surface
        nexttile
        hold on;    box  on;    grid on;    axis padded
        view(3)
        
        % surface plot
        s0 = surf(X1mesh,X2mesh,zeros(length(Z1mesh)));
            s0.LineStyle = "none";
            s0.FaceColor = 'k';
            s0.FaceAlpha = 0.2;
            
        s1 = surf(X1mesh,X2mesh,Z1mesh);
            s1.LineStyle = "none";
            s1.FaceColor = 'r';
            s1.FaceAlpha = 0.4;
        
        s2 = surf(X1mesh,X2mesh,Z2mesh);
            s2.LineStyle = "none";    
            s2.FaceColor = 'b';
            s2.FaceAlpha = 0.4;
        
        % labels
        xlabel('$x_1$');   ylabel('$x_2$');
        title('Surface plot')
        legend('$zero$','$f_1$','$f_2$')
                
        %%% contour plot of the zero level
        nexttile
        hold on;    box  on;    grid on;    axis padded
        
        % contour plot level 0
        [~,c1] = contour(X1mesh,X2mesh,Z1mesh,[0 0]);
            c1.EdgeColor = 'r';
        
        [~,c2] = contour(X1mesh,X2mesh,Z2mesh,[0 0]);
            c2.EdgeColor = 'b';
        
        % labels
        xlabel('$x_1$');   ylabel('$x_2$');  zlabel('$f(x_1,x_2)$')
        title('Contour plot')
        legend('$f_1=0$','$f_2=0$',location='northwest')
    end
    
end

%%% --- END FUNCTIONS ---
