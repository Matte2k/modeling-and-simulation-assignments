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
            nozzle.lining.mat = 'Silicon nitride';
            nozzle.lining.rho = [];             % material density          [kg/m^3]
            nozzle.lining.k   = 29;             % material conductivity     [W/(K*m)]
            nozzle.lining.l   = 0.005;          % material thickness        [m]
            nozzle.lining.c   = [];             % material specific heat    [J/(K*kg)]
            %nozzle.lining.K   = 1.724e-4;

        nozzle.conductor = struct;
            nozzle.conductor.mat = 'Tungsten';      
            nozzle.conductor.rho = 19300;       % material density          [kg/m^3]
            nozzle.conductor.k   = 164;         % material conductivity     [W/(K*m)]
            nozzle.conductor.l   = 0.030;       % material thickness        [m]
            nozzle.conductor.c   = 134;         % material specific heat    [J/(K*kg)]
            %nozzle.conductor.C   = 7.759e4; 
            %nozzle.conductor.K   = 1.829e-4;

        nozzle.interface = struct;
            nozzle.interface.mat = '';      
            nozzle.interface.rho = [];          % material density          [kg/m^3]
            nozzle.interface.k   = [];          % material conductivity     [W/(K*m)]
            nozzle.interface.l   = [];          % material thickness        [m]
            nozzle.interface.K   = 2e-3;        % material thickness        [K/W]
            nozzle.interface.c   = [];          % material specific heat    [J/(K*kg)]

        nozzle.insulator = struct;
            nozzle.insulator.mat = 'Graphene';
            nozzle.insulator.rho = 1410;        % material density          [kg/m^3]
            nozzle.insulator.k   = 0.53;        % material conductivity     [W/(K*m)]
            nozzle.insulator.l   = 0.015;       % material thickness        [m]
            nozzle.insulator.c   = 1113;        % material specific heat    [J/(K*kg)]
            %nozzle.insulator.C   = 2.354e4; 
            %nozzle.insulator.K   = 2.830e-2;

        nozzle.coating = struct;
            nozzle.coating.mat = 'Composite';
            nozzle.coating.rho = [];            % material density          [kg/m^3]
            nozzle.coating.k   = 5.2;           % material conductivity     [W/(K*m)]
            nozzle.coating.l   = 0.005;         % material thickness        [m]
            nozzle.coating.c   = [];            % material specific heat    [J/(K*kg)]
            %nozzle.coating.K   = 9.615e-4;

        nozzle.diameter = [];            % nozzle cross-section      [m^2]
        nozzle.section  = 1;             % nozzle cross-section      [m^2]
        nozzle.length   = 0.5;           % nozzle length             [m]

end
