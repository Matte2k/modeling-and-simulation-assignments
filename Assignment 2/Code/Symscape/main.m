% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment # 2
% Author: Matteo Baio 10667431

clearvars;  close all;  clc;

%% Data input

% load all system data
[time, temp, nozzle] = initData();      

% build boundary temperature function
temp.tInt = @(t) (t<1)  * (temp.tInt0 + ((temp.tInmSys.T1 - temp.tInt0).*t)/(time.mSys.T1 - time.t0)) + ...
                 (t>=1) * (temp.tInmSys.T1);
temp.tOut = @(t) temp.tExt;

% compute thermal resistance and capacitance from data
nozzle.lining    = cmpPart(nozzle.lining,    nozzle.section);
nozzle.conductor = cmpPart(nozzle.conductor, nozzle.section);
nozzle.insulator = cmpPart(nozzle.insulator, nozzle.section);
nozzle.coating   = cmpPart(nozzle.coating,   nozzle.section);

% compute mathematical model 
mModel = cmpMathModel(temp,nozzle);

%%% TASK 1
% setup a function that allows to compute and retrive results for different
% ode solver and compare with each other 
% IDEAL COMPARE: CPU time cost, step dimension (min and max), relError

%%% TASK 2
% implement if possible a way to compute eigenvalues time evolution to
% understand if there is a clear solver better than other based on the
% stability region


%% TO REFACTOR

opt = odeset('MaxStep',0.1);
[tVec,Tx] = ode89(odeSys,[time.t0,time.tf],[300,300],opt);
%[tVec,Tx] = ode113(odeSys,[time.t0,time.tf],[0,0],opt);

% ODE SOLUTION
for i = 1:length(tVec)
    Ti_ode(i,1) = Ti(tVec(i));
    To_ode(i,1) = To(tVec(i));
end

T1_ode = mModel.T1(Ti_ode,Tx(:,1));
T2_ode = Tx(:,1);
T3_ode = mModel.T3(Tx(:,1),Tx(:,2));
T4_ode = Tx(:,2);
T5_ode = mModel.T5(Tx(:,2),To_ode);

% Plot
figure(Name='input')
tiledlayout(2,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,Ti_ode)
    title('T intern')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,To_ode)
    title('T outer')
    
figure(Name='output 1')
tiledlayout(2,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,T2_ode)
    title('T 2')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,T4_ode)
    title('T 4')

figure(Name='output 2')
tiledlayout(3,1)
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,T1_ode)
    title('T 1')

    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,T3_ode)
    title('T 3')
        
    nexttile
    hold on;    grid on;    axis padded;    box on;
    plot(tVec,T5_ode)
    title('T 5')


%% Functions

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
                         1/(model.C4) * ( (1/model.Rc)*(model.T3(Tx(1),Tx(2))      - Tx(2)) ...
                                            - (1/model.Rd)*(Tx(2) - model.T5(Tx(2),model.To(t)) ) ) ];

end

