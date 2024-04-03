function [x,t,info] = bdf1(f,x0,tmax,h,bdfOptions)
%BDF1 Summary of this function goes here
%   Detailed explanation goes here
%
%   MUST ADD THE CHECK THAT: size(x0)=[n,1]

%%% Optional input definition
if nargin < 5
    bdfOptions.method = [];
    bdfOptions.alpha = [];
    bdfOptions.beta = [];
    bdfOptions.gamma = [];
end


%%% Default method
if isempty(bdfOptions.method)
    bdfOptions.method = 'Standard';
end

if isempty(bdfOptions.options)
    bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );  % default fsolve options
end


%%% Parameters definition
switch bdfOptions.method
    case 'Standard'
        alpha = 1;
        beta = 1;
        gamma = 1;

        if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                bdfOptions.method);
        end

    case 'Custom'
        if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,1]))            % TO DEBUG
            error('Parameters matrix dimensions are invalid\n');

        end

        alpha = bdfOptions.alpha;
        beta = bdfOptions.beta;
        gamma = bdfOptions.gamma;

    otherwise
        error('Insert a valid method as input\n');

end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = size(x0,1);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = [x0 zeros(dimSys,length(t)-1)];    % solution vector allocation


%%% BDF1 loop
for i = 1 : (length(t)-1)         % main loop of the method
    xk1 = x(:,i);
    
    tp = t(i) + h;
    fp = @(xp) 1/alpha * ( gamma(1)*xk1 + beta*h*f(xp,tp) ) - xp;   % k+1
    [x(:,i+1),~,conv,info] = fsolve(fp, x(:,i), bdfOptions.options);         % initial guess xk
    feval = feval + info.funcCount;

    % convergence check
    if not(isequal(conv,ones(1,length(conv))))
        warning('Implicit equation at the step %d has not been solved correctly',i)
    end
end

elapsedTime = toc(timerStart);   % timer stop

if nargout == 3
    info = struct;
        info.timeCost = elapsedTime;
        info.fevalCost = feval;
end

end