function [x,t,info] = rk1(f,x0,tmax,h,rkOptions)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

%%% Optional input definition
if nargin < 5
    rkOptions.submethod = [];
    rkOptions.alpha = [];
    rkOptions.beta = [];
end


%%% Default submethod
if isempty(rkOptions.submethod)
    rkOptions.submethod = 'FowardEuler';       % Heun parameters sets as default
end

if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))            % TO DEBUG
    warning('Parameters matrix unused, standard %s parameters are used instead\n',...
        rkOptions.submethod);
end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = length(x0);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = zeros(1,length(t));         % solution vector allocation
x(1) = x0;                      % initial condition in solution vector


%%% RK1 loop
for i=1:(length(t)-1)
    x(i+1) = x(i) + h * f(x(i),t(i));
    feval = feval + dimSys;             % function evaluation counter update
end

elapsedTime = toc(timerStart);   % timer stop

info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;

end

