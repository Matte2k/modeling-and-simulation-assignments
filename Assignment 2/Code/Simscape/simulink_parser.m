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

% build models
mModel = cmpMathModel(temp,nozzle);
sModel = cmpSimscapeModel();


function [sModel] = cmpSimscapeModel()
%
%
%
%

    % load simulink model
    sModel.name = "rocketNozzle";
    load_system(sModel.name)
       
    % set simulink solver option for current model
    %set_param(sModel.name,'Solver','auto');
    set_param(sModel.name,'StopTime','60');
    set_param(sModel.name,'MaxStep', '0.1');
    
    % execute simulink model
    sModel.simOut = sim(sModel.name);
    close_system(sModel.name,0);
    
    % retrive simulink simulation solutions for model case 1
    sModel.case1.T1 = sModel.simOut.simlog.T1_1.T.series.values('K');
    sModel.case1.T2 = sModel.simOut.simlog.T2_1.T.series.values('K');
    sModel.case1.T3 = sModel.simOut.simlog.T3_1.T.series.values('K');
    sModel.case1.T4 = sModel.simOut.simlog.T4_1.T.series.values('K');
    sModel.case1.T5 = sModel.simOut.simlog.T5_1.T.series.values('K');
    sModel.case1.time = sModel.simOut.simlog.T1_1.T.series.time;
    
    % retrive simulink simulation solutions for model case 2
    sModel.case2.T1    = sModel.simOut.simlog.T1_2.T.series.values('K');
    sModel.case2.T2in  = sModel.simOut.simlog.T2_in_2.T .series.values('K');
    sModel.case2.T2out = sModel.simOut.simlog.T2_out_2.T.series.values('K');
    sModel.case2.T3    = sModel.simOut.simlog.T3_2.T.series.values('K');
    sModel.case2.T4in  = sModel.simOut.simlog.T4_in_2.T .series.values('K');
    sModel.case2.T4out = sModel.simOut.simlog.T4_out_2.T.series.values('K');
    sModel.case2.T5    = sModel.simOut.simlog.T5_2.T.series.values('K');
    sModel.case2.time = sModel.simOut.simlog.T1_2.T.series.time;

end



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

