function [x,t,info] = rk2(f,x0,tmax,h,rkOptions)
%RK2 Summary of this function goes here
%   Detailed explanation goes here

%%% Optional input definition
if nargin < 5
    rkOptions.method = [];
    rkOptions.alpha = [];
    rkOptions.beta = [];
end


%%% Default method
if isempty(rkOptions.method)
    rkOptions.method = 'Heun';       % Heun parameters sets as default
end


%%% Parameters definition
switch rkOptions.method
    case 'Heun'
        alpha2 = [1 1]';                 % rk alpha matrix
        beta2 = [1, 0; 0.5, 0.5];        % rk beta matrix

        if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                rkOptions.method);
        end

    case 'MidPoint'
        alpha2 = [0.5 1]';                 % rk alpha matrix
        beta2 = [0.5, 0; 0, 1];            % rk beta matrix

        if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                rkOptions.method);
        end

    case 'Custom'
        if isempty(rkOptions.alpha) || isempty(rkOptions.beta)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        elseif not(isequal(size(rkOptions.alpha),[2,1])) || not(isequal(size(rkOptions.beta),[2,2]))
            error('Parameters matrix dimensions are invalid\n');
            
        end

        alpha2 = rkOptions.alpha;
        beta2 = rkOptions.beta;

        %%% Check if the custom parameters are valid
        eC1 = beta2(2,1) + beta2(2,2) == 1;
        eC2 = 2 * beta2(1,1) * beta2(2,2) == 1;
        eC3 = 2 * alpha2(1,1) * beta2(2,2) == 1;

        if eC1 ~= true ||  eC2 ~= true ||  eC3 ~= true
            error('Insert a valid parameters matrix as input');
        end

    otherwise
        error('Insert a valid method as input');
end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = size(x0,1);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = [x0 zeros(dimSys,length(t)-1)];    % solution vector allocation


%%% RK2 loop
for i=1:(length(t)-1)
    fk = f(x(:,i),t(i));
    xp = x(:,i) + beta2(1,1) * h * fk;    % x(i)=xk && t(i)=tk
    tp = t(i) + alpha2(1,1) * h;
    x(:,i+1) = x(:,i) + alpha2(2,1) * h * (beta2(2,1) * fk + beta2(2,2) * f(xp,tp));
    feval = feval + 2*dimSys;             % function evaluation counter update
end

elapsedTime = toc(timerStart);   % timer stop

info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;


end

