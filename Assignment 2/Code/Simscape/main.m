% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment # 2
% Author: Matteo Baio 10667431

clearvars;  close all;  clc;
cmap = graphicSettings();
addpath("figures\");

%% Data input

% load all system data
[time, temp, nozzle] = initData();      

% build boundary temperature function
temp.tInt = @(t) (t<1)  * (temp.tInt0 + ((temp.tInt1 - temp.tInt0).*t)/(time.t1 - time.t0)) + ...
                 (t>=1) * (temp.tInt1);
temp.tOut = @(t) temp.tExt;

% compute thermal resistance and capacitance from data
nozzle.lining    = cmpPart(nozzle.lining,    nozzle.section);
nozzle.conductor = cmpPart(nozzle.conductor, nozzle.section);
nozzle.insulator = cmpPart(nozzle.insulator, nozzle.section);
nozzle.coating   = cmpPart(nozzle.coating,   nozzle.section);

% compute mathematical model 
mModel = cmpMathModel(temp,nozzle);

odeOptions = odeset('MaxStep',0.1);
solverOpt  = solverSettings([time.t0,time.tf], [300,300], odeOptions,1);

[ode45.sol, ode45.cost]   = mModelSolver(mModel, 'ode45', solverOpt);
[ode89.sol, ode89.cost]   = mModelSolver(mModel, 'ode89', solverOpt);
[ode113.sol, ode113.cost] = mModelSolver(mModel, 'ode113', solverOpt);
[ode15s.sol, ode15s.cost] = mModelSolver(mModel, 'ode15s', solverOpt);
%[ode15i.sol, ode15i.cost] = mModelSolver(mModel, 'ode15i', solverOpt);


%%% TASK 1
% find a way to open and run Simscape simulation from matlab

%% Compare with Simscape (TEMPORARY)

%[ssData] = simscapeData(figPath,)

open('figures\T1.fig');
a = get(gca,'Children');
T1xdata = get(a,'XData');
T1ydata = get(a,'YData');
close all;

open('figures\T2.fig');
a = get(gca,'Children');
T2xdata = get(a,'XData');
T2ydata = get(a,'YData');
close all;

open('figures\T3.fig');
a = get(gca,'Children');
T3xdata = get(a,'XData');
T3ydata = get(a,'YData');
close all;

open('figures\T4.fig');
a = get(gca,'Children');
T4xdata = get(a,'XData');
T4ydata = get(a,'YData');
close all;

open('figures\T5.fig');
a = get(gca,'Children');
T5xdata = get(a,'XData');
T5ydata = get(a,'YData');
close all;


% Temperature nodes In and Out
figure(Name='input')
tiledlayout(2,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.Ti,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.Ti,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.Ti, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.Ti, 'Color', cmap(4,:))
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T intern')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.To,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.To,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.To, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.To, 'Color', cmap(4,:))
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T outer')

% Temperature nodes 2 - 4
figure(Name='output 1')
tiledlayout(2,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.T2,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.T2,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.T2, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.T2, 'Color', cmap(4,:))
    plot(T2xdata,T2ydata)
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 2')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.T4,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.T4,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.T4, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.T4, 'Color', cmap(4,:))
    plot(T4xdata,T4ydata)
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 4')

% Temperature nodes 1 - 3 - 5
figure(Name='output 2')
tiledlayout(3,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.T1,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.T1,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.T1, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.T1, 'Color', cmap(4,:))
    plot(T1xdata,T1ydata)
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 1')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.T3,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.T3,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.T3, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.T3, 'Color', cmap(4,:))
    plot(T3xdata,T3ydata)
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 3')
        
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(ode45.sol.time, ode45.sol.T5,  'Color', cmap(1,:))
    plot(ode89.sol.time, ode89.sol.T5,  'Color', cmap(2,:))
    plot(ode113.sol.time,ode113.sol.T5, 'Color', cmap(3,:))
    plot(ode15s.sol.time,ode15s.sol.T5, 'Color', cmap(4,:))
    plot(T5xdata,T5ydata)
    legend('ode45','ode89','ode113','ode15s', Location='best')
    title('T 5')


% CPU-Time cost histogram
plot.timeCostVec = [ode45.cost ode89.cost ode113.cost ode15s.cost];
plot.timeLabel   = categorical({'ode45', 'ode89','ode113','ode15s'});
plot.timeLabel   = reordercats(plot.timeLabel,{'ode45', 'ode89','ode113','ode15s'});
figure(Name='time cost')
    hold on;    grid on;    axis padded;    box on;
    bar(plot.timeLabel,plot.timeCostVec)
    ylabel('Time cost [s]');


%% Eigenanalysis

lambda = eig(mModel.eigSys);
    xlambda = real(lambda);
    ylambda = imag(lambda);


%% FUNCTIONS

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


function [nodesTemp, timeCost, rawOdeSol] = mModelSolver(mModel, solverName, solverOpt)
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
            timeCost = toc/solverOpt.iter;

        case 'ode113'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode113(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            timeCost = toc/solverOpt.iter;

        case 'ode78'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode78(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            timeCost = toc/solverOpt.iter;

        case 'ode89'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode89(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            timeCost = toc/solverOpt.iter;

        case 'ode15s'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode15s(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            timeCost = toc/solverOpt.iter;

        case 'ode15i'
            for i = 1:solverOpt.iter
                [rawOdeSol.tVec,rawOdeSol.Tx] = ode15i(mModel.odeSys,solverOpt.time,solverOpt.x0,solverOpt.odeOpt);
            end
            timeCost = toc/solverOpt.iter;

        otherwise
            error('select a valide ode solver');
    end


    % Nodes temperature solutions
    nodesTemp.time = rawOdeSol.tVec;
    for i = 1:length(rawOdeSol.tVec)
        nodesTemp.Ti(i,1) = mModel.Ti(nodesTemp.time(i));
        nodesTemp.To(i,1) = mModel.To(nodesTemp.time(i));
    end

    nodesTemp.T1 = mModel.T1(nodesTemp.Ti,rawOdeSol.Tx(:,1));
    nodesTemp.T2 = rawOdeSol.Tx(:,1);
    nodesTemp.T3 = mModel.T3(rawOdeSol.Tx(:,1),rawOdeSol.Tx(:,2));
    nodesTemp.T4 = rawOdeSol.Tx(:,2);
    nodesTemp.T5 = mModel.T5(rawOdeSol.Tx(:,2),nodesTemp.To);


end

