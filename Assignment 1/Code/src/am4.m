function [x,t,info] = am4(f,x0,tmax,h,amOptions)
%AM4 Summary of this function goes here
%   Detailed explanation goes here



%%% Optional input definition
if nargin < 5
    amOptions.method = [];
    amOptions.alpha = [];
    amOptions.beta = [];
end


%%% Default method
if isempty(amOptions.method)
    amOptions.method = 'Standard';       % Heun parameters sets as default
end

if isempty(amOptions.options)
    amOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );  % default fsolve options
end


%%% Parameters definition
switch amOptions.method
    case 'Standard'
        alpha = 24;
        beta = [9 19 -5 1];

        if not(isempty(amOptions.alpha)) || not(isempty(amOptions.beta))            % TO DEBUG
            warning('Parameters matrix unused, standard %s parameters are used instead\n',...
                amOptions.method);
        end

    case 'Custom'
        if isempty(amOptions.alpha) || isempty(amOptions.beta)            % TO DEBUG
            error('No parameters matrix has been given as input\n');

        elseif not(isequal(size(amOptions.alpha),[1,1])) || not(isequal(size(amOptions.beta),[1,4]))            % TO DEBUG
            error('Parameters matrix dimensions are invalid\n');

        end

        alpha = amOptions.alpha;
        beta = amOptions.beta;

    otherwise
        error('Insert a valid method as input\n');

end


%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = size(x0,1);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
x = [x0 zeros(dimSys,length(t)-4)];    % solution vector allocation


%%% AM1 loop
for i = 4 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    
    tp = t(i) + h;
    fp = @(xp) x(:,i) + h/alpha * (beta(1)*f(xp,tp) + beta(2)*fk1 + beta(3)*fk2 + beta(4)*fk3) - xp;   % k+1
    [x(:,i+1),~,conv,info] = fsolve(fp, x(:,i), amOptions.options);         % initial guess xk
    feval = feval + info.funcCount + 3*dimSys; 

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