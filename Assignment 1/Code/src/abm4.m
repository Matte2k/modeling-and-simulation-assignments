function [x,t,info] = abm4(f,x0,tmax,h,abmOptions)
%ABM4 Summary of this function goes here
%   Detailed explanation goes here



%%% Optional input definition
if nargin < 5
    abmOptions.method = [];
    abmOptions.alpha = [];
    abmOptions.betaP = [];
    abmOptions.betaC = [];
end


%%% Default method
if isempty(abmOptions.method)
    abmOptions.method = 'Standard';       % Heun parameters sets as default
end

%%% Parameters definition
switch abmOptions.method
    case 'Standard'
        alpha = 24;
        betaP = [55 -59  37  -9];
        betaC = [9   19  -5   1];

        if not(isempty(abmOptions.alpha)) || not(isempty(abmOptions.betaP)) || not(isempty(abmOptions.betaC))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                abmOptions.method);
        end

    case 'Custom'
        if isempty(abmOptions.alpha) || isempty(abmOptions.betaP) || isempty(abmOptions.betaC)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        elseif not(isequal(size(abmOptions.alpha),[1,1])) || not(isequal(size(abmOptions.betaP),[1,4])) || not(isequal(size(abmOptions.betaC),[1,4]))            % TO DEBUG
            error('Parameters matrix dimensions are invalid\n');

        end

        alpha = abmOptions.alpha;
        betaP = abmOptions.betaP;
        betaC = abmOptions.betaC;

    otherwise
        error('Insert a valid method as input\n');

end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = size(x0,1);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = [x0 zeros(dimSys,length(t)-4)];    % solution vector allocation


%%% ABM4 loop
for i = 4 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    fk4 = f(x(:,i-3),t(i-3));
    xp = x(:,i) + h/alpha * (betaP(1)*fk1 + betaP(2)*fk2 + betaP(3)*fk3 + betaP(4)*fk4);

    tp = t(i) + h;
    x(:,i+1) = x(:,i) + h/alpha * (betaC(1)*f(xp,tp) + betaC(2)*fk1 + betaC(3)*fk2 + betaC(4)*fk3);   % k+1
    feval = feval + 5*dimSys; 
end

elapsedTime = toc(timerStart);   % timer stop

if nargout == 3
info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;
end

end