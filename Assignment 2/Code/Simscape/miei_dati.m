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
        time.vec = linspace(time.t0,time.tf,500);

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