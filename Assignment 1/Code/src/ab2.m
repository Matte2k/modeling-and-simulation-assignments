function [x,t,info] = ab2(f,x0,tmax,h,abOptions)
%AB2 Summary of this function goes here
%   Detailed explanation goes here

%%% Optional input definition
if nargin < 5
    abOptions.submethod = [];
    abOptions.alpha = [];
    abOptions.beta = [];
end


%%% Default submethod
if isempty(abOptions.submethod)
    abOptions.submethod = 'Standard';       % Heun parameters sets as default
end

%%% Parameters definition
switch abOptions.submethod
    case 'Standard'
        alpha = 2;
        beta = [3 -1];

        if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                abOptions.submethod);
        end

    case 'Custom'
        if isempty(abOptions.alpha) || isempty(abOptions.beta)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

            % elseif size(rkOptions.alpha) ~= size(zeros(2,1)) || size(rkOptions.beta) ~= size(zeros(2,2))
            %     error('Parameters matrix dimensions are invalid\n');   
            % ISSUE:  the 'size' operation don't give a scalar as result so
            %         there are problems in logical operations
        end

        alpha = abOptions.alpha;
        beta = abOptions.beta;

    otherwise
        error('Insert a valid submethod as input\n');

end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = length(x0);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = zeros(dimSys,length(t));         % solution vector allocation
x(:,1) = x0;                      % initial condition in solution vector


%%% AB2 loop
for i = 2 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    x(:,i+1) = x(:,i) + h/alpha(2) * (beta(2,1)*fk1 + beta(2,2)*fk2);
    feval = feval + 2*dimSys;             % function evaluation counter update
end

elapsedTime = toc(timerStart);   % timer stop

info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;

end

