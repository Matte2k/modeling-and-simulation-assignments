function [x,t,info] = ab3(f,x0,tmax,h,abOptions)
%AB3 Summary of this function goes here
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
        alpha = 12;
        beta = [23 -16  5];

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


%%% AB3 loop
for i = 2 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    x(:,i+1) = x(:,i) + h/alpha(3) * (beta(3,1)*fk1 + beta(3,2)*fk2 + beta(3,3)*fk3);
    feval = feval + 3*dimSys;             % function evaluation counter update
end

elapsedTime = toc(timerStart);   % timer stop

info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;

end
