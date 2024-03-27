function [x,t,info] = rk4(f,x0,tmax,h,rkOptions)
%RK4 Summary of this function goes herefIVP
%   Detailed explanation goes here

%%% Optional input definition
if nargin < 5
    rkOptions.submethod = [];
    rkOptions.alpha = [];
    rkOptions.beta = [];
end


%%% Default submethod
if isempty(rkOptions.submethod)
    rkOptions.submethod = 'Runge-Kutta';       % Heun parameters sets as default
end


%%% Parameters definition
switch rkOptions.submethod
    case 'Runge-Kutta'
        alpha4 = [0.5 0.5 1 1]';           % rk alpha matrix
        beta4 = diag([0.5 0.5 1 1/6]);     % rk beta matrix
        beta4(4,:) = [1/6 1/3 1/3 1/6];

        if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                rkOptions.submethod);
        end

    case 'Butcher'                  % NOT IMPLEMENTED YET
        %alpha4 = []';           % TBD
        %beta4 = diag([]);     % rk beta matrix
        %beta4(4,:) = [];

        if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                rkOptions.submethod);
        end
        error ('Method has not been implemented yet\n');

    case 'Custom'
        if isempty(rkOptions.alpha) || isempty(rkOptions.beta)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        % elseif size(rkOptions.alpha) ~= size(zeros(2,1)) || size(rkOptions.beta) ~= size(zeros(2,2))
        %     error('Parameters matrix dimensions are invalid\n');
        end

        alpha4 = rkOptions.alpha;
        beta4 = rkOptions.beta;

    otherwise
        error('Insert a valid submethod as input\n');
end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = length(x0);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = zeros(dimSys,length(t));    % solution vector allocation
x(:,1) = x0;                      % initial condition in solution vector


%%% RK4 loop
for i=1:(length(t)-1)                                % calculation loop
    fk = f(x(:,i),t(i));
    xp1 = x(:,i) + beta4(1,1) * h * fk;    % x(i)=xk && t(i)=tk
    tp1 = t(i) + alpha4(1,1) * h;

    xp2 = x(:,i) + beta4(2,2) * h * f(xp1,tp1);    % x(i)=xk && t(i)=tk
    tp2 = t(i) + alpha4(2,1) * h;

    xp3 = x(:,i) + beta4(3,3) * h * f(xp2,tp2);    % x(i)=xk && t(i)=tk
    tp3 = t(i) + alpha4(3,1) * h;

    x(:,i+1) = x(:,i) + alpha4(4,1) * h * ( beta4(4,1) * fk ...
        + beta4(4,2) * f(xp1,tp1) ...
        + beta4(4,3) * f(xp2,tp2)...
        + beta4(4,4) * f(xp3,tp3) );
    
    feval = feval + 4*dimSys;             % function evaluation counter update
end

elapsedTime = toc(timerStart);   % timer stop

info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;

end