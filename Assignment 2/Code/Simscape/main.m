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
% compute thermal resistance and capacitance from data
nozzle.lining    = cmpPart(nozzle.lining,    nozzle.section);
nozzle.conductor = cmpPart(nozzle.conductor, nozzle.section);
nozzle.insulator = cmpPart(nozzle.insulator, nozzle.section);
nozzle.coating   = cmpPart(nozzle.coating,   nozzle.section);

% compute mathematical model 
mModel = cmpMathModel(temp,nozzle);

odeOptions = odeset('MaxStep',0.1);
solverOpt  = solverSettings([time.t0,time.tf], [300,300], odeOptions, CPUtimeCostIter);

mModel.ode45  = mModelSolver(mModel, 'ode45', solverOpt);
mModel.ode89  = mModelSolver(mModel, 'ode89', solverOpt);
mModel.ode113 = mModelSolver(mModel, 'ode113', solverOpt);
mModel.ode15s = mModelSolver(mModel, 'ode15s', solverOpt);
%mModel.ode15i = mModelSolver(mModel, 'ode15i', solverOpt);

% compute simscape model
sModel.auto = execSimscapeModel('rocketNozzle', CPUtimeCostIter, caseModelSelector, [time.t0,time.tf]);


%% Result Plot

% Temperature nodes In and Out
figure(Name='input')
tiledlayout(2,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(mModel.ode45.time,  mModel.ode45.Ti,  'Color', cmap(1,:))
    plot(mModel.ode89.time,  mModel.ode89.Ti,  'Color', cmap(2,:))
    plot(mModel.ode113.time, mModel.ode113.Ti, 'Color', cmap(3,:))
    plot(mModel.ode15s.time, mModel.ode15s.Ti, 'Color', cmap(4,:))
    plot(sModel.auto.case1.time, sModel.auto.case1.Ti , 'Color', cmap(5,:))
    legend('ode45','ode89','ode113','ode15s','Simscape', Location='best')
    title('T intern')

    nexttile
    hold on;    grid on;    axis padded;    box on;

    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T outer')

% Temperature nodes 2 - 4
figure(Name='output 1')
tiledlayout(2,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;


    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 2')

    nexttile
    hold on;    grid on;    axis padded;    box on;


    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 4')

% Temperature nodes 1 - 3 - 5
figure(Name='output 2')
tiledlayout(3,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;

    plot(sModel.auto.case1.time, sModel.auto.case1.T1 , 'Color', cmap(5,:))
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 1')

    nexttile
    hold on;    grid on;    axis padded;    box on;


    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 3')

    nexttile
    hold on;    grid on;    axis padded;    box on;


    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 5')


% CPU-Time cost histogram
visual.timeCostVec = [mModel.ode45.CPUtime mModel.ode89.CPUtime mModel.ode113.CPUtime mModel.ode15s.CPUtime];
visual.timeLabel   = categorical({'ode45', 'ode89','ode113','ode15s'});
visual.timeLabel   = reordercats(visual.timeLabel,{'ode45', 'ode89','ode113','ode15s'});
figure(Name='time cost')
    hold on;    grid on;    axis padded;    box on;
    bar(visual.timeLabel,visual.timeCostVec)
    ylabel('Time cost [s]');


% Eigenvalues
figure(Name='Eigenvalues')
    hold on;    grid on;    axis equal;    box on;
    plot(real(mModel.lambda), imag(mModel.lambda) , 'o', MarkerFaceColor='auto')
    yline(0)
    xline(0)



%% FUNCTIONS

function [cmap] = graphicSettings()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Graphic settings
    set(0,'defaulttextinterpreter','latex');  
    set(0,'defaultAxesTickLabelInterpreter','latex');  
    set(0,'defaultLegendInterpreter','latex');
    set(0,'DefaultLineLineWidth',1.2)

    cmap = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

    warning('off','MATLAB:fplot:NotVectorized');
    
end


function [time, temp, nozzle] = initData()
%INITIALIZATION DATA - initialize all the data for the simulation
%
%   Syntax:     
%
%
%   Reference:  MatWeb
%

    time = struct;
        time.t0 = 0;        % initial simulation time       [sec]
        time.t1 = 1;        % end of fire-up time           [sec]
        time.tf = 60;       % final simulation time         [sec]
    
    temp = struct;
        temp.tInt0 =  293.15;       % initial internal temperature:     20 C°  |   293.15 K
        temp.tInt1 = 1273.15;       % regime  internal temperature:   1000 C°  |  1273.15 K
        temp.tExt  =  293.15;       % fixed   external tempertaure:     20 C°  |   293.15 K
    
    nozzle = struct;   
        nozzle.lining = struct;
            nozzle.lining.mat = 'Inconel';
            nozzle.lining.rho = 8.440;          % material density          [kg/m^3]
            nozzle.lining.k   = 14.5;           % material conductivity     [W/(K*m)]
            nozzle.lining.l   = 0.32;           % material thickness        [m]
            nozzle.lining.c   = [];             % material specific heat    [J/(K*kg)]
    
        nozzle.conductor = struct;
            nozzle.conductor.mat = 'Copper';      
            nozzle.conductor.rho = 8.960;       % material density          [kg/m^3]
            nozzle.conductor.k   = 400;         % material conductivity     [W/(K*m)]
            nozzle.conductor.l   = 0.54;        % material thickness        [m]
            nozzle.conductor.c   = 385;         % material specific heat    [J/(K*kg)]
    
        nozzle.interface = struct;
            nozzle.interface.mat = '';      
            nozzle.interface.rho = [];          % material density          [kg/m^3]
            nozzle.interface.k   = [];          % material conductivity     [W/(K*m)]
            nozzle.interface.l   = 0.18;        % material thickness        [m]
            nozzle.interface.K   = 0.000678;    % material thickness        [K/W]
            nozzle.interface.c   = [];          % material specific heat    [J/(K*kg)]
    
        nozzle.insulator = struct;
            nozzle.insulator.mat = 'Ceramic';
            nozzle.insulator.rho = 2.800;       % material density          [kg/m^3]
            nozzle.insulator.k   = 1.5;         % material conductivity     [W/(K*m)]
            nozzle.insulator.l   = 0.36;        % material thickness        [m]
            nozzle.insulator.c   = 800;         % material specific heat    [J/(K*kg)]
    
        nozzle.coating = struct;
            nozzle.coating.mat = 'Stainless Steel';
            nozzle.coating.rho = 8.000;         % material density          [kg/m^3]
            nozzle.coating.k   = 16;            % material conductivity     [W/(K*m)]
            nozzle.coating.l   = 0.36;          % material thickness        [m]
            nozzle.coating.c   = [];            % material specific heat    [J/(K*kg)]
        
        nozzle.section = 1.3273;                % nozzle cross-section      [m^2]
        nozzle.length  = nozzle.lining.l ...    % nozzle length             [m]
                          + nozzle.conductor.l ...
                          + nozzle.interface.l ... 
                          + nozzle.insulator.l ...
                          + nozzle.coating.l;               

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
    model.R1 = (nozzle.lining.K);
    model.Ra = (nozzle.lining.K    + nozzle.conductor.K) * 1/2;
    model.Rb = (nozzle.conductor.K + nozzle.interface.K) * 1/2;
    model.Rc = (nozzle.interface.K + nozzle.insulator.K) * 1/2;
    model.Rd = (nozzle.insulator.K + nozzle.coating.K) * 1/2;
    model.R5 = (nozzle.coating.K);
    
    % Mathematical system thermal capacitance
    model.C2 = nozzle.conductor.C;
    model.C4 = nozzle.insulator.C;
    
    % Boundary temperature
    model.Ti = @(t) temp.tInt(t);
    model.To = @(t) temp.tOut(t);
    
    % Nodes temperature function
    model.T1 = @(Ti,T2) Ti*((2*model.Ra)/(2*model.Ra+model.R1)) + T2*(model.R1/(2*model.Ra+model.R1));
    model.T3 = @(T2,T4) T2*(model.Rc/(model.Rc+model.Rb)) + T4*(model.Rb/(model.Rc+model.Rb));
    model.T5 = @(T4,To) T4*(model.R5/(model.R5+2*model.Rd)) + To*((2*model.Rd)/(model.R5+2*model.Rd));
    
    % Differential system   [mSys.Tx(1)=T2  and  mSys.Tx(2)=T4]
    model.odeSys = @(t,Tx) [ 1/(model.C2) * ( (1/model.Ra)*(model.T1(model.Ti(t),Tx(1)) - Tx(1)) ...
                                            - (1/model.Rb)*(Tx(1) - model.T3(Tx(1),Tx(2)     ) ) )  ; ...
                             1/(model.C4) * ( (1/model.Rc)*(model.T3(Tx(1),Tx(2))       - Tx(2)) ...
                                            - (1/model.Rd)*(Tx(2) - model.T5(Tx(2),model.To(t)) ) ) ];
    % model.odeSys = @(t,Tx) [ 1/(model.C2) * ( Tx(1)*(model.R1/(model.Ra*(2*model.Ra + model.R1)) - 1/model.Ra - 1/model.Rb + model.Rc/((model.Rc + model.Rb)*model.Rb) ) ...
    %                                         + Tx(2)*(1/(model.Rc+model.Rb)) + model.Ti(t)*(2/(2*model.Ra+model.R1))  )  ; ...
    %                          1/(model.C4) * ( Tx(2)*(model.Rb/(model.Rc*(model.Rc + model.Rb)) - 1/model.Rc - 1/model.Rd + model.R5/(model.Rd*(model.R5 + 2*model.Rd)) ) ...
    %                                         + Tx(1)*(1/(model.Rc+model.Rb)) + model.To(t)*(2/(2*model.Rd+model.R5))  ) ];

    % State space matrix of the differential system (for eigenanalysis)
    model.eigSys = [ 1/(model.C2) * (model.R1/(model.Ra*(2*model.Ra + model.R1)) - 1/model.Ra - 1/model.Rb + model.Rc/((model.Rc + model.Rb)*model.Rb) ), ...
                     1/(model.C2) * (1/(model.Rc+model.Rb)) ; ...
                     1/(model.C4) * (model.Rb/(model.Rc*(model.Rc + model.Rb)) - 1/model.Rc - 1/model.Rd + model.R5/(model.Rd*(model.R5 + 2*model.Rd)) ), ...
                     1/(model.C4) * (1/(model.Rc+model.Rb)) ];
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
        set_param(modelName,'StartTime', '0.1');
    else
        set_param(modelName,'StopTime',  num2str(timeLimit(2)));
        set_param(modelName,'StartTime', num2str(timeLimit(1)));
    end

    if nargin >= 5 && ~isempty(solverName)
        set_param(modelName,'Solver', solverName);
    end

    if nargin >= 6 && ~isempty(maxStep)
        set_param(modelName,'MaxStep', maxStep);
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
        modelSol.case2.To = simulinkModel.simOut.simlog.To_2.T.series.values('K');
        modelSol.case2.T1    = simulinkModel.simOut.simlog.T1_2.T.series.values('K');
        modelSol.case2.T1    = simulinkModel.simOut.simlog.T1_2.T.series.values('K');
        modelSol.case2.T2in  = simulinkModel.simOut.simlog.T2_in_2.T .series.values('K');
        modelSol.case2.T2out = simulinkModel.simOut.simlog.T2_out_2.T.series.values('K');
        modelSol.case2.T3    = simulinkModel.simOut.simlog.T3_2.T.series.values('K');
        modelSol.case2.T4in  = simulinkModel.simOut.simlog.T4_in_2.T .series.values('K');
        modelSol.case2.T4out = simulinkModel.simOut.simlog.T4_out_2.T.series.values('K');
        modelSol.case2.T5    = simulinkModel.simOut.simlog.T5_2.T.series.values('K');
        modelSol.case2.To    = simulinkModel.simOut.simlog.To_2.A.T.series.values('K');
        modelSol.case2.time = simulinkModel.simOut.simlog.T1_2.T.series.time;
    end

    % export all the simulation log
    modelSol.simulationOutput = simulinkModel.simOut;
    modelSol.modelName = modelName;

    % close simulink model
    close_system(modelName,0);

end
