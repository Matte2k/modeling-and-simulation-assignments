function [x,t,info] = bdf4(f,x0,tmax,h,bdfOptions)
%BDF4 Summary of this function goes here
%   Detailed explanation goes here
%
%   MUST ADD THE CHECK THAT: size(x0)=[n,4]


%%% Optional input definition
if nargin < 5
    bdfOptions.method = [];
    bdfOptions.alpha = [];
    bdfOptions.beta = [];
    bdfOptions.gamma = [];
end


%%% Default method
if isempty(bdfOptions.method)
    bdfOptions.method = 'Standard';       % Heun parameters sets as default
end

if isempty(bdfOptions.options)
    bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );  % default fsolve options
end


%%% Parameters definition
switch bdfOptions.method
    case 'Standard'
        alpha = 25;
        beta = 12;
        gamma = [48 -36 16 3];

        if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                bdfOptions.method);
        end

    case 'Custom'
        if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,4]))            % TO DEBUG
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
x = [x0 zeros(dimSys,length(t)-4)];    % solution vector allocation


%%% AB4 loop
for i = 4 : (length(t)-1)         % main loop of the method
    xk1 = x(:,i  );
    xk2 = x(:,i-1);
    xk3 = x(:,i-2);
    xk4 = x(:,i-3);

    tp = t(i) + h;
    fp = @(xp) 1/alpha * ( gamma(1)*xk1 + gamma(2)*xk2 + gamma(3)*xk3 + gamma(4)*xk4 + beta*h*f(xp,tp) ) - xp;   % k+1
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