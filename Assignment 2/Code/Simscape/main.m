% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment # 2
% Author: Matteo Baio 10667431

clearvars;  close all;  clc;
cmap = graphicSettings();


%% Data input

% load all system data
CPUtimeCostIter   = 1;
caseModelSelector = 1;
[time, temp, nozzle] = initData(); 

% build boundary temperature function
temp.tInt = @(t) (t<1)  * (temp.tInt0 + ((temp.tInt1 - temp.tInt0).*t)/(time.t1 - time.t0)) + ...
                 (t>=1) * (temp.tInt1);
temp.tOut = @(t) temp.tExt;


%% Model build and simulation
%compute thermal resistance and capacitance from data
nozzle.lining    = cmpPart(nozzle.lining,    nozzle.section);
nozzle.conductor = cmpPart(nozzle.conductor, nozzle.section);
nozzle.insulator = cmpPart(nozzle.insulator, nozzle.section);
nozzle.coating   = cmpPart(nozzle.coating,   nozzle.section);

% compute mathematical model 
mModel = cmpMathModel(temp,nozzle);

odeOptions = odeset('MaxStep',0.01);
solverOpt  = solverSettings([time.t0,time.tf], [293.15 293.15], odeOptions , CPUtimeCostIter);

mModel.ode45  = mModelSolver(mModel, 'ode45',  solverOpt);
mModel.ode89  = mModelSolver(mModel, 'ode89',  solverOpt);
mModel.ode113 = mModelSolver(mModel, 'ode113', solverOpt);
mModel.ode23t = mModelSolver(mModel, 'ode23t', solverOpt);
mModel.ode15s = mModelSolver(mModel, 'ode15s', solverOpt);

% ---------------------------------------------
% new ode solver
F = ode;
F.ODEFcn = mModel.odeSys;
F.InitialValue = [293.15 293.15];
F.AbsoluteTolerance = 1e-12;
F.RelativeTolerance = 1e-6;
F.Solver = "auto";
sol = solve(F,0,60);
F.SelectedSolver
% ---------------------------------------------

% compute simscape model
sModel.auto = execSimscapeModel('rocketNozzle', CPUtimeCostIter, caseModelSelector, [time.t0,time.tf],'ode23t',0.01);

% interpolate simscape and ode data on same grid
mModel.ode45.interpT1 = interp1(mModel.ode45.time, mModel.ode45.T1, time.vec);
mModel.ode45.interpT2 = interp1(mModel.ode45.time, mModel.ode45.T2, time.vec);
mModel.ode45.interpT3 = interp1(mModel.ode45.time, mModel.ode45.T3, time.vec);
mModel.ode45.interpT4 = interp1(mModel.ode45.time, mModel.ode45.T4, time.vec);
mModel.ode45.interpT5 = interp1(mModel.ode45.time, mModel.ode45.T5, time.vec);

mModel.ode89.interpT1 = interp1(mModel.ode89.time, mModel.ode89.T1, time.vec);
mModel.ode89.interpT2 = interp1(mModel.ode89.time, mModel.ode89.T2, time.vec);
mModel.ode89.interpT3 = interp1(mModel.ode89.time, mModel.ode89.T3, time.vec);
mModel.ode89.interpT4 = interp1(mModel.ode89.time, mModel.ode89.T4, time.vec);
mModel.ode89.interpT5 = interp1(mModel.ode89.time, mModel.ode89.T5, time.vec);

mModel.ode23t.interpT1 = interp1(mModel.ode23t.time, mModel.ode23t.T1, time.vec);
mModel.ode23t.interpT2 = interp1(mModel.ode23t.time, mModel.ode23t.T2, time.vec);
mModel.ode23t.interpT3 = interp1(mModel.ode23t.time, mModel.ode23t.T3, time.vec);
mModel.ode23t.interpT4 = interp1(mModel.ode23t.time, mModel.ode23t.T4, time.vec);
mModel.ode23t.interpT5 = interp1(mModel.ode23t.time, mModel.ode23t.T5, time.vec);

mModel.ode15s.interpT1 = interp1(mModel.ode15s.time, mModel.ode15s.T1, time.vec);
mModel.ode15s.interpT2 = interp1(mModel.ode15s.time, mModel.ode15s.T2, time.vec);
mModel.ode15s.interpT3 = interp1(mModel.ode15s.time, mModel.ode15s.T3, time.vec);
mModel.ode15s.interpT4 = interp1(mModel.ode15s.time, mModel.ode15s.T4, time.vec);
mModel.ode15s.interpT5 = interp1(mModel.ode15s.time, mModel.ode15s.T5, time.vec);

mModel.ode113.interpT1 = interp1(mModel.ode113.time, mModel.ode113.T1, time.vec);
mModel.ode113.interpT2 = interp1(mModel.ode113.time, mModel.ode113.T2, time.vec);
mModel.ode113.interpT3 = interp1(mModel.ode113.time, mModel.ode113.T3, time.vec);
mModel.ode113.interpT4 = interp1(mModel.ode113.time, mModel.ode113.T4, time.vec);
mModel.ode113.interpT5 = interp1(mModel.ode113.time, mModel.ode113.T5, time.vec);

sModel.auto.case1.interpT1 = interp1(sModel.auto.case1.time, sModel.auto.case1.T1, time.vec);
sModel.auto.case1.interpT2 = interp1(sModel.auto.case1.time, sModel.auto.case1.T2, time.vec);
sModel.auto.case1.interpT3 = interp1(sModel.auto.case1.time, sModel.auto.case1.T3, time.vec);
sModel.auto.case1.interpT4 = interp1(sModel.auto.case1.time, sModel.auto.case1.T4, time.vec);
sModel.auto.case1.interpT5 = interp1(sModel.auto.case1.time, sModel.auto.case1.T5, time.vec);

% compute error between simscape and ode solver
error.ode45.T1 = sModel.auto.case1.interpT1 - mModel.ode45.interpT1;
error.ode45.T2 = sModel.auto.case1.interpT2 - mModel.ode45.interpT2;
error.ode45.T3 = sModel.auto.case1.interpT3 - mModel.ode45.interpT3;
error.ode45.T4 = sModel.auto.case1.interpT4 - mModel.ode45.interpT4;
error.ode45.T5 = sModel.auto.case1.interpT5 - mModel.ode45.interpT5;

error.ode89.T1 = sModel.auto.case1.interpT1 - mModel.ode89.interpT1;
error.ode89.T2 = sModel.auto.case1.interpT2 - mModel.ode89.interpT2;
error.ode89.T3 = sModel.auto.case1.interpT3 - mModel.ode89.interpT3;
error.ode89.T4 = sModel.auto.case1.interpT4 - mModel.ode89.interpT4;
error.ode89.T5 = sModel.auto.case1.interpT5 - mModel.ode89.interpT5;

error.ode23t.T1 = sModel.auto.case1.interpT1 - mModel.ode23t.interpT1;
error.ode23t.T2 = sModel.auto.case1.interpT2 - mModel.ode23t.interpT2;
error.ode23t.T3 = sModel.auto.case1.interpT3 - mModel.ode23t.interpT3;
error.ode23t.T4 = sModel.auto.case1.interpT4 - mModel.ode23t.interpT4;
error.ode23t.T5 = sModel.auto.case1.interpT5 - mModel.ode23t.interpT5;

error.ode15s.T1 = sModel.auto.case1.interpT1 - mModel.ode15s.interpT1;
error.ode15s.T2 = sModel.auto.case1.interpT2 - mModel.ode15s.interpT2;
error.ode15s.T3 = sModel.auto.case1.interpT3 - mModel.ode15s.interpT3;
error.ode15s.T4 = sModel.auto.case1.interpT4 - mModel.ode15s.interpT4;
error.ode15s.T5 = sModel.auto.case1.interpT5 - mModel.ode15s.interpT5;

error.ode113.T1 = sModel.auto.case1.interpT1 - mModel.ode113.interpT1;
error.ode113.T2 = sModel.auto.case1.interpT2 - mModel.ode113.interpT2;
error.ode113.T3 = sModel.auto.case1.interpT3 - mModel.ode113.interpT3;
error.ode113.T4 = sModel.auto.case1.interpT4 - mModel.ode113.interpT4;
error.ode113.T5 = sModel.auto.case1.interpT5 - mModel.ode113.interpT5;


%% Result Plot

% Temperature nodes In and Out
figure(Name='Input nodes')
tiledlayout(1,2)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.Ti,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.Ti,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.Ti, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.Ti, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.Ti , 'Color', cmap(5,:))
    %legend('ode45','ode89','ode113','ode15s','Simscape', Location='best')
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]');
    title('Temperature inner lining')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.To,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.To,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.To, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.To, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.To , 'Color', cmap(5,:))
    %legend('ode45','ode89','ode113','ode15s','Simscape', Location='best')
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]');
    title('Temperature outer lining')
set(gcf,'units','centimeters','position',[0,0,20,8]);
%exportgraphics(gcf,'figure\ex1_tempInput.eps');
  

% Temperature nodes 2 - 4
figure(Name='Unknown nodes')
tiledlayout(1,2)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.T2,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.T2,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.T2, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.T2, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T2 , 'Color', cmap(5,:))
    %legend('ode45','ode89','ode113','ode15s','Simscape', Location='best')
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Temperature node 2')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.T4,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.T4,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.T4, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.T4, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T4 , 'Color', cmap(5,:))
    %legend('ode45','ode89','ode113','ode15s','Simscape', Location='best')
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Temperature node 4')
set(gcf,'units','centimeters','position',[0,0,20,8]);
%exportgraphics(gcf,'figure\ex1_tempKeyNodes.eps');

% Temperature nodes 1 - 3 - 5
figure(Name='Remaining nodes')
tiledlayout(1,3)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.T1,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.T1,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.T1, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.T1, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T1 , 'Color', cmap(5,:))
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Temperature node 1')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.T3,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.T3,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.T3, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.T3, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T3 , 'Color', cmap(5,:))
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Temperature node 3')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.T5,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.T5,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.T5, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.T5, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T5 , 'Color', cmap(5,:))
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Temperature node 5')

lgd = legend('ode45','ode89','ode113','ode15s','Simscape');
lgd.Layout.Tile = 'south';
lgd.Orientation = 'horizontal';
set(gcf,'units','centimeters','position',[0,0,20,8]);
%exportgraphics(gcf,'figure\ex1_tempOtherNodes.eps');

% Overall temperature compare case 1
figure(Name='all temp case 1')
tiledlayout(1,2)
nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode23t.time, mModel.ode23t.T1 , 'Color', cmap(1,:))
    plot(mModel.ode23t.time, mModel.ode23t.T2 , 'Color', cmap(2,:))
    plot(mModel.ode23t.time, mModel.ode23t.T3 , 'Color', cmap(3,:))
    plot(mModel.ode23t.time, mModel.ode23t.T4 , 'Color', cmap(4,:))
    plot(mModel.ode23t.time, mModel.ode23t.T5 , 'Color', cmap(5,:))
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Matlab simulation')

nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(sModel.auto.case1.time, sModel.auto.case1.T1 , 'Color', cmap(1,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T2 , 'Color', cmap(2,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T3 , 'Color', cmap(3,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T4 , 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.T5 , 'Color', cmap(5,:))
    xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
    title('Simscape simulation')
lgd = legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$');
lgd.Layout.Tile = 'south';
lgd.Orientation = 'horizontal';    
set(gcf,'units','centimeters','position',[0,0,20,9]);
%exportgraphics(gcf,'figure\ex1_tempCompare.eps');

% Overall temperature compare case 2
% figure(Name='all temp case 2')
%     hold on;    grid on;    axis padded;    box on;
%     plot(sModel.auto.case2.time, sModel.auto.case2.T1 , 'Color', cmap(1,:))
%     plot(sModel.auto.case2.time, sModel.auto.case2.T2 , 'Color', cmap(2,:))
%     plot(sModel.auto.case2.time, sModel.auto.case2.T3 , 'Color', cmap(3,:))
%     plot(sModel.auto.case2.time, sModel.auto.case2.T4 , 'Color', cmap(4,:))
%     plot(sModel.auto.case2.time, sModel.auto.case2.T5 , 'Color', cmap(5,:))
%     legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$', Location='best')
%     xlabel('Time [$sec$]');   ylabel('Temperature [$K$]')
%     title('Simscape simulation')
% set(gcf,'units','centimeters','position',[0,0,11,9]);
% %exportgraphics(gcf,'figure\ex1_tempCase2.eps');

% CPU-Time cost histogram
visual.timeCostVec = [mModel.ode45.CPUtime mModel.ode89.CPUtime mModel.ode113.CPUtime mModel.ode23t.CPUtime mModel.ode15s.CPUtime sModel.auto.CPUtime];
visual.timeLabel   = categorical({'ode45', 'ode89','ode113','ode23t','ode15s','simscape'});
visual.timeLabel   = reordercats(visual.timeLabel,{'ode45','ode89','ode113','ode23t','ode15s','simscape'});
figure(Name='time cost')
    hold on;    grid on;    axis padded;    box on;
    bar(visual.timeLabel,visual.timeCostVec)
    ylabel('Time cost [s]');    xlabel('Integration method');
set(gcf,'units','centimeters','position',[0,0,11,9]);
%exportgraphics(gcf,'figure\ex1_CPUtime.eps');

% Eigenvalues
figure(Name='Eigenvalues')
    hold on;    grid on;    axis equal;    box on;
    yline(0,'LineWidth',1);  ylim([-2e-3 , 2e-3]);
    xline(0,'LineWidth',1);  xlim([-10e-3, 1e-3]);
    plot(real(mModel.lambda(1)), imag(mModel.lambda(1)), 'o', MarkerFaceColor='r',MarkerEdgeColor='r')
    plot(real(mModel.lambda(2)), imag(mModel.lambda(2)), 'o', MarkerFaceColor='b',MarkerEdgeColor='b')
    xlabel('$Re_{h\lambda}$');     ylabel('$Im_{h\lambda}$');
    legend('','','$\lambda_1$','$\lambda_2$');
set(gcf,'units','centimeters','position',[0,0,11,6]);
%exportgraphics(gcf,'figure/ex1_eigenvalues.eps');

% Error
figure(Name='Error ode89')
    hold on;    grid on;    axis padded;    box on;
    plot(time.vec,abs(error.ode89.T1))
    plot(time.vec,abs(error.ode89.T2))
    plot(time.vec,abs(error.ode89.T3))
    plot(time.vec,abs(error.ode89.T4))
    plot(time.vec,abs(error.ode89.T5))
    legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$', Location='best')

figure(Name='Error ode23t')
    hold on;    grid on;    axis padded;    box on;
    plot(time.vec,abs(error.ode23t.T1))
    plot(time.vec,abs(error.ode23t.T2))
    plot(time.vec,abs(error.ode23t.T3))
    plot(time.vec,abs(error.ode23t.T4))
    plot(time.vec,abs(error.ode23t.T5))
    legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$', Location='best')

figure(Name='Error ode45')
    hold on;    grid on;    axis padded;    box on;
    plot(time.vec,abs(error.ode45.T1))
    plot(time.vec,abs(error.ode45.T2))
    plot(time.vec,abs(error.ode45.T3))
    plot(time.vec,abs(error.ode45.T4))
    plot(time.vec,abs(error.ode45.T5))
    legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$', Location='best')

figure(Name='Error ode15s')
    hold on;    grid on;    axis padded;    box on;
    plot(time.vec,abs(error.ode15s.T1))
    plot(time.vec,abs(error.ode15s.T2))
    plot(time.vec,abs(error.ode15s.T3))
    plot(time.vec,abs(error.ode15s.T4))
    plot(time.vec,abs(error.ode15s.T5))
    legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$', Location='best')


%% FUNCTIONS

function [cmap] = graphicSettings()
%GRAPHIC SETTINGS - set graphics options for the visual output
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


% function [time, temp, nozzle] = initData()      % quelli del Bo
% %INITIALIZATION DATA - initialize all the data for the simulation
% %
% %   Syntax:     
% %
% %
% %   Reference:  MatWeb
% %
% 
%     time = struct;
%         time.t0  = 0;        % initial simulation time       [sec]
%         time.t1  = 1;        % end of fire-up time           [sec]
%         time.tf  = 60;       % final simulation time         [sec]
%         time.vec = time.t0:1:time.tf;
% 
%     temp = struct;
%         temp.tInt0 =  293.15;       % initial internal temperature:     20 C°  |   293.15 K
%         temp.tInt1 = 1273.15;       % regime  internal temperature:   1000 C°  |  1273.15 K
%         temp.tExt  =  293.15;       % fixed   external tempertaure:     20 C°  |   293.15 K
% 
%     nozzle = struct;   
%         nozzle.lining = struct;
%             nozzle.lining.mat = 'Silicon nitride';
%             nozzle.lining.rho = [];             % material density          [kg/m^3]
%             nozzle.lining.k   = 29;             % material conductivity     [W/(K*m)]
%             nozzle.lining.l   = 0.005;          % material thickness        [m]
%             nozzle.lining.c   = [];             % material specific heat    [J/(K*kg)]
%             %nozzle.lining.K   = 1.724e-4;
% 
%         nozzle.conductor = struct;
%             nozzle.conductor.mat = 'Tungsten';      
%             nozzle.conductor.rho = 19300;       % material density          [kg/m^3]
%             nozzle.conductor.k   = 164;         % material conductivity     [W/(K*m)]
%             nozzle.conductor.l   = 0.030;       % material thickness        [m]
%             nozzle.conductor.c   = 134;         % material specific heat    [J/(K*kg)]
%             %nozzle.conductor.C   = 7.759e4; 
%             %nozzle.conductor.K   = 1.829e-4;
% 
%         nozzle.interface = struct;
%             nozzle.interface.mat = '';      
%             nozzle.interface.rho = [];          % material density          [kg/m^3]
%             nozzle.interface.k   = [];          % material conductivity     [W/(K*m)]
%             nozzle.interface.l   = [];          % material thickness        [m]
%             nozzle.interface.K   = 2e-3;        % material thickness        [K/W]
%             nozzle.interface.c   = [];          % material specific heat    [J/(K*kg)]
% 
%         nozzle.insulator = struct;
%             nozzle.insulator.mat = 'Graphene';
%             nozzle.insulator.rho = 1410;        % material density          [kg/m^3]
%             nozzle.insulator.k   = 0.53;        % material conductivity     [W/(K*m)]
%             nozzle.insulator.l   = 0.015;       % material thickness        [m]
%             nozzle.insulator.c   = 1113;        % material specific heat    [J/(K*kg)]
%             %nozzle.insulator.C   = 2.354e4; 
%             %nozzle.insulator.K   = 2.830e-2;
% 
%         nozzle.coating = struct;
%             nozzle.coating.mat = 'Composite';
%             nozzle.coating.rho = [];            % material density          [kg/m^3]
%             nozzle.coating.k   = 5.2;           % material conductivity     [W/(K*m)]
%             nozzle.coating.l   = 0.005;         % material thickness        [m]
%             nozzle.coating.c   = [];            % material specific heat    [J/(K*kg)]
%             %nozzle.coating.K   = 9.615e-4;
% 
%         nozzle.diameter = [];            % nozzle cross-section      [m^2]
%         nozzle.section  = 1;             % nozzle cross-section      [m^2]
%         nozzle.length   = 0.5;           % nozzle length             [m]
% 
% end

function [time, temp, nozzle] = initData()
%INITIALIZATION DATA - initialize all the data for the simulation
%
%   Syntax:     
%
%
%   Reference:  MatWeb
%

    time = struct;
        time.t0  = 0;        % initial simulation time       [sec]
        time.t1  = 1;        % end of fire-up time           [sec]
        time.tf  = 60;       % final simulation time         [sec]
        time.vec = time.t0:1:time.tf;

    temp = struct;
        temp.tInt0 =  293.15;       % initial internal temperature:     20 C°  |   293.15 K
        temp.tInt1 = 1273.15;       % regime  internal temperature:   1000 C°  |  1273.15 K
        temp.tExt  =  293.15;       % fixed   external tempertaure:     20 C°  |   293.15 K

    nozzle = struct;   
        nozzle.lining = struct;
            nozzle.lining.mat = 'Inconel';
            nozzle.lining.rho = 8440;           % material density          [kg/m^3]
            nozzle.lining.k   = 14.5;           % material conductivity     [W/(K*m)]
            nozzle.lining.l   = 0.02;           % material thickness        [m]
            nozzle.lining.c   = [];             % material specific heat    [J/(K*kg)]

        nozzle.conductor = struct;
            nozzle.conductor.mat = 'Copper';      
            nozzle.conductor.rho = 8960;        % material density          [kg/m^3]
            nozzle.conductor.k   = 400;         % material conductivity     [W/(K*m)]
            nozzle.conductor.l   = 0.03;        % material thickness        [m]
            nozzle.conductor.c   = 385;         % material specific heat    [J/(K*kg)]

        nozzle.interface = struct;
            nozzle.interface.mat = '';      
            nozzle.interface.rho = [];          % material density          [kg/m^3]
            nozzle.interface.k   = [];          % material conductivity     [W/(K*m)]
            nozzle.interface.l   = 0.01;        % material thickness        [m]
            nozzle.interface.K   = 0.000678;    % material thickness        [K/W]
            nozzle.interface.c   = [];          % material specific heat    [J/(K*kg)]

        nozzle.insulator = struct;
            nozzle.insulator.mat = 'Ceramic';
            nozzle.insulator.rho = 2800;        % material density          [kg/m^3]
            nozzle.insulator.k   = 1.5;         % material conductivity     [W/(K*m)]
            nozzle.insulator.l   = 0.02;        % material thickness        [m]
            nozzle.insulator.c   = 800;         % material specific heat    [J/(K*kg)]

        nozzle.coating = struct;
            nozzle.coating.mat = 'Stainless Steel';
            nozzle.coating.rho = 8000;          % material density          [kg/m^3]
            nozzle.coating.k   = 16;            % material conductivity     [W/(K*m)]
            nozzle.coating.l   = 0.02;          % material thickness        [m]
            nozzle.coating.c   = [];            % material specific heat    [J/(K*kg)]

        nozzle.section  = 1;                    % nozzle cross-section      [m^2]
        nozzle.length   = 0.5;                  % nozzle length             [m]
              

end

function [part] = cmpPart(part,crossSection)
%COMPUTE PART - compute thermal resistence and/or capacitance of the part
%if possible
%
%   COMPLETE DESCRIPTION HERE
%
%

    % Compute thermal resistence of the part if possible
    if ~isempty(part.k)
        part.K = part.l / (part.k * crossSection);
    end
    
    % Compute thermal capacity of the part if possible
    if ~isempty(part.c)
        part.C = part.c * (part.rho * part.l * crossSection);
    end

end


function [model] = cmpMathModel(temp,nozzle)
%COMPUTE MATHEMATICAL MODEL - compute the mathematical model of the system
%
%   COMPLETE DESCRIPTION HERE
%
%

    % Mathematical system thermal resistence
    model.Ri2 = (nozzle.lining.K) + (nozzle.conductor.K * 1/2) ;
    model.R24 = (nozzle.conductor.K * 1/2) + nozzle.interface.K + (nozzle.insulator.K * 1/2);
    model.R4o = (nozzle.insulator.K * 1/2) + (nozzle.coating.K);
    
    % Mathematical system thermal capacitance
    model.C2 = nozzle.conductor.C;
    model.C4 = nozzle.insulator.C;
    
    % Mathematical system heat flow
    model.Qi2 = @(Ti,T2) (Ti-T2)/model.Ri2;
    model.Q24 = @(T2,T4) (T2-T4)/model.R24;
    model.Q4o = @(T4,To) (T4-To)/model.R4o;

    % Boundary temperature
    model.Ti = @(t) temp.tInt(t);
    model.To = @(t) temp.tOut(t);
    
    % Differential system   [mSys.Tx(1)=T2  and  mSys.Tx(2)=T4]
    model.odeSys = @(t,Tx) [ 1/(model.C2) * ( (model.Ti(t)-Tx(1))/model.Ri2 - (Tx(1)-Tx(2))/model.R24 ) ; ...
                             1/(model.C4) * ( (Tx(1)-Tx(2))/model.R24 - (Tx(2)-model.To(t))/model.R4o ) ];
    
    % Nodes temperature function
    model.T1 = @(Ti,T2) Ti - model.Qi2(Ti,T2) * (nozzle.lining.K/2);
    model.T3 = @(T2,T4) T4 + model.Q24(T2,T4) * ((nozzle.interface.K + nozzle.insulator.K) * 1/2);
    model.T5 = @(T4,To) To + model.Q4o(T4,To) * (nozzle.coating.K/2);

    % State space matrix of the differential system (for eigenanalysis)
    model.eigSys = [ -1/(model.C2*model.Ri2)-1/(model.C2*model.R24) , ...
                      1/(model.C2*model.R24) ; ...
                      1/(model.C4*model.R24) , ...
                     -1/(model.C4*model.R24)-1/(model.C4*model.R4o) ];
    model.lambda = eig(model.eigSys);

end


function [solverOpt] = solverSettings(timeLimit, initGuess, odeOptions, odeIter)
%
%   COMPLETE DESCRIPTION HERE
%
%
    if nargin < 4 || isempty(odeIter)
        odeIter = 1;
    end

    if length(timeLimit) == 1
        timeLimit(2) = timeLimit;
        timeLimit(1) = 0;
    end

    solverOpt = struct;
        solverOpt.time   = [timeLimit(1),timeLimit(2)];
        solverOpt.x0     = initGuess;
        solverOpt.odeOpt = odeOptions;
        solverOpt.iter   = odeIter;

end


function [modelSol, rawOdeSol] = mModelSolver(mModel, solverName, solverOpt)
%
%
%
%

    switch solverName
        case 'ode45'
            tic
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode45(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        case 'ode113'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode113(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        case 'ode78'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode78(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        case 'ode89'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode89(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        case 'ode15s'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode15s(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        case 'ode15i'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode15i(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        case 'ode23t'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode23t(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            modelSol.CPUtime = toc/solverOpt.iter;

        otherwise
            error('select a valide ode solver');
    end

    % Nodes temperature solutions
    modelSol.time = rawOdeSol.tVec;
    for i = 1:length(rawOdeSol.tVec)
        modelSol.Ti(i,1) = mModel.Ti(modelSol.time(i));
        modelSol.To(i,1) = mModel.To(modelSol.time(i));
    end

    modelSol.T1 = mModel.T1(modelSol.Ti,rawOdeSol.Tx(:,1));
    modelSol.T2 = rawOdeSol.Tx(:,1);
    modelSol.T3 = mModel.T3(rawOdeSol.Tx(:,1),rawOdeSol.Tx(:,2));
    modelSol.T4 = rawOdeSol.Tx(:,2);
    modelSol.T5 = mModel.T5(rawOdeSol.Tx(:,2),modelSol.To);

end


function [modelSol] = execSimscapeModel(modelName, timeCostIter, caseID, timeLimit, solverName, maxStep)
%
%
%
%
    
    % case model selector
    if nargin < 3 || isempty(caseID)
        caseID = 0;
    end

    % load simulink model
    load_system(modelName)
       
    % set simulink solver option for current model
    if nargin < 4 || isempty(timeLimit)
        set_param(modelName,'StopTime',  '60');
        set_param(modelName,'StartTime', '0');
    else
        set_param(modelName,'StopTime',  num2str(timeLimit(2)));
        set_param(modelName,'StartTime', num2str(timeLimit(1)));
    end

    if nargin >= 5 && ~isempty(solverName)
        set_param(modelName,'Solver', solverName);
    end

    if nargin >= 6 && ~isempty(maxStep)
        set_param(modelName,'MaxStep', num2str(maxStep));
    end

    % execute simulink model
    tic
    for i = 1:timeCostIter
        simulinkModel.simOut = sim(modelName);
    end
    modelSol.CPUtime = toc/timeCostIter;
    
    if caseID == 1 || caseID == 0
        % retrive simulink simulation solutions for model case 1
        modelSol.case1.Ti = simulinkModel.simOut.simlog.Ti_1.T.series.values('K');
        modelSol.case1.T1 = simulinkModel.simOut.simlog.T1_1.T.series.values('K');
        modelSol.case1.T2 = simulinkModel.simOut.simlog.T2_1.T.series.values('K');
        modelSol.case1.T3 = simulinkModel.simOut.simlog.T3_1.T.series.values('K');
        modelSol.case1.T4 = simulinkModel.simOut.simlog.T4_1.T.series.values('K');
        modelSol.case1.T5 = simulinkModel.simOut.simlog.T5_1.T.series.values('K');
        modelSol.case1.To = simulinkModel.simOut.simlog.To_1.A.T.series.values('K');
        modelSol.case1.time = simulinkModel.simOut.simlog.T1_1.T.series.time;
    end

    if caseID == 2 || caseID == 0
        % retrive simulink simulation solutions for model case 2
        modelSol.case2.To    = simulinkModel.simOut.simlog.To_2.T.series.values('K');
        modelSol.case2.T1    = simulinkModel.simOut.simlog.T1_2.T.series.values('K');
        modelSol.case2.T1    = simulinkModel.simOut.simlog.T1_2.T.series.values('K');
        modelSol.case2.T2in  = simulinkModel.simOut.simlog.T2_in_2.T .series.values('K');
        modelSol.case2.T2out = simulinkModel.simOut.simlog.T2_out_2.T.series.values('K');
        modelSol.case2.T3    = simulinkModel.simOut.simlog.T3_2.T.series.values('K');
        modelSol.case2.T4in  = simulinkModel.simOut.simlog.T4_in_2.T .series.values('K');
        modelSol.case2.T4out = simulinkModel.simOut.simlog.T4_out_2.T.series.values('K');
        modelSol.case2.T5    = simulinkModel.simOut.simlog.T5_2.T.series.values('K');
        modelSol.case2.To    = simulinkModel.simOut.simlog.To_2.A.T.series.values('K');
        modelSol.case2.time  = simulinkModel.simOut.simlog.T1_2.T.series.time;
    end

    % export all the simulation log
    modelSol.simulationOutput = simulinkModel.simOut;
    modelSol.modelName = modelName;

    % close simulink model
    close_system(modelName,0);

end


