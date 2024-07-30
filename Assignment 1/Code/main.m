% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment # 1
% Author: Matteo Baio 10667431

addpath(genpath('src/'));
addpath(genpath('figure/'));

%% Ex 1
%code ex1
clearvars; close all; clc;
graphicSettings;

%%% DATA INPUT
f = @(x1,x2) [x2.^2 - x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
xGuess = [[3.5 2.5]',[2.5 -2]'];
toll = 1e-12;
nmax = 1e3;
config = struct;
    config.print = true;
    config.plot = false;


%%% INITIAL GUESS DEFINITION
% Initial guess graphic analysis
f1 = @(x1,x2) x2.^2 - x1 - 2;           f2 = @(x1,x2) -x1.^2 + x2 + 10;
[Z1mesh,Z2mesh] = zerosGuess(f1,f2);            % plot functions
plot(xGuess(1,:),xGuess(2,:),'ko');             % plot guess
legend('$f_1=0$','$f_2=0$','Initial guess')     % plot legend
set(gcf,'units','centimeters','position',[0,0,20,10]);
exportgraphics(gcf,'figure\ex1_initGuess.eps',Resolution=1000);


% Variable initialization
numZeros = size(xGuess,2);
solutionSym = zeros(2,numZeros);    solutionFD  = zeros(2,numZeros);    solutionCD  = zeros(2,numZeros);
convergeSym = zeros(1,numZeros);    convergeFD =  zeros(1,numZeros);    convergeCD =  zeros(1,numZeros);
infoSym = cell(1,numZeros);         infoFD = cell(1,numZeros);          infoCD = cell(1,numZeros);


%%% NEWTON'S METHODS
for i = 1:numZeros
    [solutionSym(:,i),convergeSym(i),infoSym{i}] = newton(f, xGuess(:,i), 's', toll, nmax, config); % symbolic math
    [solutionFD(:,i), convergeFD(i), infoFD{i}]  = newton(f, xGuess(:,i), 'f', toll, nmax, config); % forward difference
    [solutionCD(:,i), convergeCD(i), infoCD{i}]  = newton(f, xGuess(:,i), 'c', toll, nmax, config); % centered difference

    % Error compare between different submethod
    errSym = infoSym{i}.errorVector(end);
    errFD  = infoFD{i}.errorVector(end);
    errCD  = infoCD{i}.errorVector(end);
    
    % Compare error to find best method
    if errSym < errFD && errSym < errCD
        fprintf('Sym is the best method for the zero %d \n\n',i);
    elseif errFD < errSym && errFD < errCD
        fprintf('FD is the best method for the zero %d \n\n',i);
    elseif errCD < errFD && errCD < errSym
        fprintf('CD is the best method for the zero %d \n\n',i);
    end
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
        xlabel('$Iteration$');    ylabel('$Error$');
        legend('$FD$','$CD$','$sym$',Location='best')
        tString = sprintf('$z_{%d}=[%.4f,%.4f]$',i,solutionFD(i,1),solutionFD(i,1));
        title(tString)
end
set(gcf,'units','centimeters','position',[0,0,20,10]);
exportgraphics(gcf,'figure\ex1_convergence.eps');



%% Ex2 
%code ex2
clearvars; close all; clc;
[cmap]=graphicSettings();

%%% DATA INPUT
f = @(x,t)(x - 2.*(t).^2 + 2);
solutionIVP = @(t)(2.*t.^2 + 4.*t - exp(t) + 2);
avgTimeIter = 100;
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
% -add graphic export here-


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
% -add graphic export here-


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

title('Integration solution');      legend('RK2','RK4','Location','best');  % title + legend
xlabel('Final error');   ylabel('Time cost');                               % axis label

drawnow
% -add graphic export here-


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
    xlim([-4 1])
    
leg = legend('$RK2$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
% -add graphic export here-

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
    xlim([-4 1])

leg = legend('$RK4$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
% -add graphic export here-

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
% -add graphic export here-


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
hVec  = zeros(length(orders),length(tollVec));
feval = zeros(length(orders),length(tollVec));
xF    = zeros(1,length(alphaVec));
yF    = zeros(1,length(alphaVec));


%%% COMPUTE SOLUTION RK1,2,4 
for ord = 1:length(RKcell) 
    RK = RKcell{ord};
    
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
    % -add graphic export here-
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
    % -add graphic export here-


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
        xlim([-5.5 10.5])
    
    lstring{i} = sprintf('$BI2_{%.1f}$',thetaVec(i));
    leg = legend(lstring{i},'$Stable$','Orientation', 'Horizontal');
        leg.Layout.Tile = 'south';        leg.Color = 'w';
    % -add graphic export here-
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
% -add graphic export here-
    

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
    ylim([-0.2 1.2])
    
nexttile
    axis padded;     grid on;    box on;     hold on
    plot(RK4.t,RK4.x(2,:),'Color',cmap(1,:))
    plot(IEX4.t,IEX4.x(2,:),'Color',cmap(2,:))
    %plot(RK4mod.t,RK4mod.x(2,:),'--','Color',cmap(3,:))
    plot(Analytic.t,Analytic.x(2,:),'o','Color',cmap(5,:),MarkerSize=2)
    xlabel('$Time$');   ylabel('$x_2$');
    ylim([-0.2 1.2])

leg = legend('$RK4$','$IEX4$','$Analytical$','Orientation', 'Horizontal');
%leg = legend('$RK4$','$IEX4$','$RK4mod$','$Analytical$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
% -add graphic export here-
    

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
% -add graphic export here-

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
% -add graphic export here-

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
    xlim([-45 15]);
    legend('$IEX4$','$RK4$','$h\lambda$','$h_{mod}\lambda$')
    title('$RK4$ vs $IEX4$')
% -add graphic export here-


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


%%% EIGENVALUE ENVELOP
timeVec = 0:h:tmax;
figure("Name",'Eigenvalue evolution')
    hold on;    grid on;    box on;
    plot(timeVec,h*eig1(timeVec),'r')
    plot(timeVec,h*eig2(timeVec),'b')
    limitAB3 = yline(-0.6,Color=cmap(3,:),LineWidth=1.5,LineStyle='--',Label='AB3');
        limitAB3.LabelHorizontalAlignment='center';
    limitAM3 = yline(-6,Color=cmap(5,:),LineWidth=1.5,LineStyle='--',Label='AM3');
        limitAM3.LabelHorizontalAlignment='center';
    limitABM3 = yline(-1.7,Color=cmap(6,:),LineWidth=1.5,LineStyle='--',Label='ABM3');
        limitABM3.LabelHorizontalAlignment='center';        
    xlabel('Time');     ylabel('$h\lambda_{i}$');     ylim([-6.5,0.5])
    legend('$h\lambda_{1}$','$h\lambda_{2}$',Location='best')

%%% --- END CODE --- 

%% FUNCTIONS
