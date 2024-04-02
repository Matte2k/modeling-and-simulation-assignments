function [x,t,info] = ab1(f,x0,tmax,h,abOptions)
%AB1 Summary of this function goes here
%   Detailed explanation goes here
%
%   MUST ADD THE CHECK THAT: size(x0)=[n,1]

%%% Optional input definition
if nargin < 5
    abOptions.method = [];
    abOptions.alpha = [];
    abOptions.beta = [];
end


%%% Default method
if isempty(abOptions.method)
    abOptions.method = 'Standard';
end


%%% Parameters definition
switch abOptions.method
    case 'Standard'
        alpha = 1;
        beta = 1;

        if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                abOptions.method);
        end

    case 'Custom'
        if isempty(abOptions.alpha) || isempty(abOptions.beta)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        elseif not(isequal(size(abOptions.alpha),[1,1])) || not(isequal(size(abOptions.beta),[1,1]))            % TO DEBUG
            error('Parameters matrix dimensions are invalid\n');

        end

        alpha = abOptions.alpha;
        beta = abOptions.beta;

    otherwise
        error('Insert a valid method as input\n');

end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = size(x0,1);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = [x0 zeros(dimSys,length(t)-1)];    % solution vector allocation


%%% AB1 loop
for i = 1 : (length(t)-1)         % main loop of the method
    fk = f(x(:,i),t(i));
    x(:,i+1) = x(:,i) + h/alpha * (beta*fk);
    feval = feval + dimSys;             % function evaluation counter update
end

elapsedTime = toc(timerStart);   % timer stop

if nargout == 3
    info = struct;
        info.timeCost = elapsedTime;
        info.fevalCost = feval;
end

end