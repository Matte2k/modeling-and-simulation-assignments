% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment # 1
% Author: Matteo Baio 10667431

addpath(genpath('src\'));

%% Ex 1
%code ex1
clearvars; close all; clc;

%%% DATA INPUT
f = @(x1,x2) [x2.^2 - x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
toll = 1e-12;
nmax = 1e3;
config = struct;
    config.print = true;
    config.plot = true;


%%% INITIAL GUESS DEFINITION
% Initial guess graphic analysis
f1 = @(x1,x2) x2.^2 - x1 - 2;
f2 = @(x1,x2) -x1.^2 + x2 + 10;
[Z1mesh,Z2mesh] = zerosGuess(f1,f2);
% -add graphic export here-


% Variable initialization
xGuess = [[3.5 2.5]',[2.5 -2]'];
numZeros = size(xGuess,2);
solutionSym = zeros(2,numZeros);    solutionFD  = zeros(2,numZeros);    solutionCD  = zeros(2,numZeros);
convergeSym = zeros(1,numZeros);    convergeFD =  zeros(1,numZeros);    convergeCD =  zeros(1,numZeros);
infoSym = cell(1,numZeros);         infoFD = cell(1,numZeros);          infoCD = cell(1,numZeros);


%%% NEWTON'S METHODS
figure('Name','Convergence');   % Figure initialization
    t1 = tiledlayout(1,2);
    title(t1,'Newton convergence')
for i = 1:numZeros
    nexttile
        % Solution with different submethod
        [solutionSym(:,i),convergeSym(i),infoSym{i}] = newton(f, xGuess(:,i), 's', toll, nmax, config);
        hold on;
        [solutionFD(:,i), convergeFD(i), infoFD{i}]  = newton(f, xGuess(:,i), 'f', toll, nmax, config);
        [solutionCD(:,i), convergeCD(i), infoCD{i}]  = newton(f, xGuess(:,i), 'c', toll, nmax, config);
        grid on;    axis padded;    box on
        legend('$sym$','$FD$','$CD$',Location='best')
        tString = sprintf('$z_{%d}=[%.4f,%.4f]$',i,solutionFD(1,i),solutionFD(1,i));
        title(tString)
    % -add graphic export here-
    
    % Error compare between different submethod
    errSym = infoSym{i}.errorVector(end);
    errFD  = infoFD{i}.errorVector(end);
    errCD  = infoCD{i}.errorVector(end);
    
    if errSym < errFD && errSym < errCD
        fprintf('Sym is the best method for the zero %d \n\n',i);
    elseif errFD < errSym && errFD < errCD
        fprintf('FD is the best method for the zero %d \n\n',i);
    elseif errCD < errFD && errCD < errSym
        fprintf('CD is the best method for the zero %d \n\n',i);
    end

end

%%% Plot della progressione?

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
avgTimeIter = 10;
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
        semilogy(RK2{i}.t,RK2{i}.error)
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
    title(t2,'RK4 method')

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
        semilogy(RK4{i}.t,RK4{i}.error)
        hold on
    end
    grid on;    box on;    axis padded;
    title('Integration error')
    legend(lString{1:length(h)},'Location','best')
    xlabel('$Time$');   ylabel('$Error$');
 
drawnow
% -add graphic export here-


%%% POST PROCESSING   //WIP                 NEED HELP FOR THE LEGEND
% Error vs CPU time
colorPalette = ['r','g','b','m'];           % Temporary
figure('Name','Trade off')
hold on;    grid on;    box on;     axis padded
for i = 1:length(h)
    scatter(RK2{i}.error(end),RK2{i}.info.avgTimeCost,[],colorPalette(i),"filled","o");
    scatter(RK4{i}.error(end),RK4{i}.info.avgTimeCost,[],colorPalette(i),"filled","square")
    %scatterString = {scatterString, lString{i}};
end
title('Integration solution')
legend(lString{1:end},'Location','best')    % Find better way to do
xlabel('Final error');   ylabel('Time cost');

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

dimSys = 2;
stepNum = 100;

%%% RK2 STABILITY PROBLEM
figure('Name','RK2');
    t1 = tiledlayout(1,2);
    title(t1,'RK2 stability analysis')

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
    title(t2,'RK4 stability analysis')

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

%%% INPUT
x0      = [1 1]';
tmax    = 1;
tollVec = [1e-3 1e-4 1e-5 1e-6];
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
stabEdge = cell(length(orders),length(tollVec));
hVec  = zeros(length(orders),length(tollVec));
feval = zeros(length(orders),length(tollVec));
xF    = zeros(1,length(alphaVec));
yF    = zeros(1,length(alphaVec));


%%% STABILITY PLOT RK1,2,4 
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
            [hVec(ord,tols),~,conv]=fzero(G,guess);                    % Problem solution

            if i == 1
                Gcell{ord,tols} = G;                                % Initial guess analysis
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
        
        stabEdge{ord,tols} = [xF,flip(xF,2);yF,-flip(yF,2)];    % Stability region
    end
end


%%% STABILITY PLOT RK1,2,4
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
        title(t,tstring);
    
    ax1 = nexttile;
        hold (ax1,'on')
        grid on;    axis padded;    box on;
        yline(0, LineStyle = '--')
        xlim(view.ax1.x(ord,:));    ylim(view.ax1.y(ord,:));
        xlabel('$h$');              ylabel('$G(h)$');
        title('Stability function in alpha = pi')

    ax2 = nexttile;
        hold (ax2,'on')
        grid on;    axis padded;    box on;
        xline(0,'--');              yline(0,'--')
        xlim(view.ax2.x(ord,:));    ylim(view.ax2.y(ord,:));
        xlabel('$Re\{h\lambda\}$'); ylabel('$Im\{h\lambda\}$');
        title('Stability region')
    
    cmap = colormap(hsv(length(tollVec)));     % colormap definition
    
    for tols = 1:length(tollVec)
        fplot(ax1,Gcell{ord,tols},[0 hMax],'Color',cmap(tols,:));
        xline(ax1,hVec(ord,tols),'Color',cmap(tols,:),'LineStyle',':');
        fill (ax2,stabEdge{ord,tols}(1,:), stabEdge{ord,tols}(2,:),cmap(tols,:),'FaceAlpha',0.7);
    end
    
    ax1.Layer = 'top';      hold(ax1,'off');
    ax2.Layer = 'top';      hold(ax2,'off');
    
    cbar = colorbar;
        cbar.Ticks = linspace(0.1,0.9,length(tollVec));
        cbar.TickLabels = tollVec;
        cbar.Label.String = '$toll$';
        cbar.Label.Interpreter = 'latex';
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
        title(t,tstring);
    
    ax1 = nexttile;
        hold (ax1,'on')
        grid on;    axis padded;    box on;
        yline(0, LineStyle = '--')
        xlim(view.ax1.xZoom(ord,:));    ylim(view.ax1.yZoom(ord,:));
        xlabel('$h$');                  ylabel('$G(h)$');
        title('Stability function in alpha = pi')

    ax2 = nexttile;
        hold (ax2,'on')
        grid on;    axis padded;    box on;
        xline(0,'--');                  yline(0,'--')
        xlim(view.ax2.xZoom(ord,:));    ylim(view.ax2.yZoom(ord,:));
        xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
        title ('Stability region')
    
    cmap = colormap(hsv(length(tollVec)));     % colormap definition
    
    for tols = 1:length(tollVec)
        fplot(ax1,Gcell{ord,tols},[0 hMax],'Color',cmap(tols,:));
        xline(ax1,hVec(ord,tols),'Color',cmap(tols,:),'LineStyle',':');
        fill (ax2,stabEdge{ord,tols}(1,:), stabEdge{ord,tols}(2,:),cmap(tols,:),'FaceAlpha',0.7);
    end
ax1.Layer = 'top';      hold(ax1,'off');
ax2.Layer = 'top';      hold(ax2,'off');
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

thetaVec = [0.1 0.3 0.4 0.7 0.9];
dimSys   = 2;
stepNum  = 100;
stabEdge = cell(1,length(thetaVec));
lstring  = cell(1,length(thetaVec));
pal = 'grmcb';

%%% BI2 STABILITY PROBLEM
for i = 1:length(thetaVec)
    figure('Name','BI2');
        t = tiledlayout(1,2);
        tstring = sprintf('BI2_{%.1f} stability analysis',thetaVec(i));
        title(t,tstring)
    
    nexttile
        [F,guess,degVec,hVec]=selectOp('BI2',dimSys,thetaVec(i));
        stabGuess(F,degVec(1),hVec);
    
    nexttile
        [~,~,~,stabEdge{i}]=stabRegion(F,stepNum,degVec,guess);
        xlim([-5.5 10.5])
    
    lstring{i} = sprintf('$BI2_{%.1f}$',thetaVec(i));
    leg = legend(lstring{i},'$Stable$','Orientation', 'Horizontal');
        leg.Layout.Tile = 'south';        leg.Color = 'w';
    % -add graphic export here-
end

figure('Name','BI2 compare');
    axis equal;     grid on;    box on;     hold on
    cmapBI2 = colormap(hsv(length(thetaVec)));     % colormap definition
        
    for i = 1:length(thetaVec)
        plot(stabEdge{i}(1,:),stabEdge{i}(2,:),'Color',cmapBI2(i,:))
    end

cbarBI2 = colorbar;
    cbarBI2.Ticks = linspace(0.1,0.9,length(thetaVec));
    cbarBI2.TickLabels = thetaVec;
    cbarBI2.Label.String = '$h$';
    cbarBI2.Label.Interpreter = 'latex';

xline(0,LineStyle='--');        yline(0,LineStyle='--');
xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
xlim([-6 11]);
tstring = sprintf('BI2_{n} stability analysis');
title(tstring)
% -add graphic export here-
    

%% Ex6
%code ex6
clearvars; close all; clc;

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
for i = 2 : (length(Analytic.t))       % CHIEDERE A RIC SE É IL MODO PIU' SMART
    tNow = Analytic.t(i);
    Analytic.x(:,i) = solutionIVP(tNow);
end


%%% INTEGRATION RESULTS                 [Fix linestyle and color]
figure('Name','Integration results');
    t1 = tiledlayout(1,2);
    title(t1,'Integration results')

nexttile
    axis padded;     grid on;    box on;     hold on
    plot(RK4.t,RK4.x(1,:))
    plot(RK4mod.t,RK4mod.x(1,:))
    plot(IEX4.t,IEX4.x(1,:))
    plot(Analytic.t,Analytic.x(1,:),'r--')
    xlabel('$Time$');   ylabel('$x1$');
    ylim([-0.5 1.5])
    
nexttile
    axis padded;     grid on;    box on;     hold on
    plot(RK4.t,RK4.x(2,:))
    plot(RK4mod.t,RK4mod.x(2,:))
    plot(IEX4.t,IEX4.x(2,:))
    plot(Analytic.t,Analytic.x(2,:),'r--')
    xlabel('$Time$');   ylabel('$x2$');
    ylim([-0.5 1.5])

leg = legend('$RK4$','$RK4mod$','$IEX4$','$Analytical$','Orientation', 'Horizontal');
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
    title(t1,'IEX4 stability analysis')

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
    title(t2,'RK4 stability analysis')

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
    axis equal;     grid on;    box on;     hold on
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
    
    xline(0,LineStyle='--');        yline(0,LineStyle='--');
    xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
    xlim([-45 15]);
    legend('$IEX4$','$RK4$','$h\lambda\$','$h_{mod}\lambda\$')
    title('$RK4$ vs $IEX4$')
% -add graphic export here-


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

