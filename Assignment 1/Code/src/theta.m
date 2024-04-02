function [x,t,info] = theta(f,x0,tmax,h,tOptions,visualConfig)
%THETA Summary of this function goes here
%   Detailed explanation goes here

%%% Optional input definition
if nargin < 6
    visualConfig.plot = true;
end

if nargin < 5
tOptions = struct;
    tOptions.theta = [];
    tOptions.method = [];
    tOptions.submethod = [];
        tOptions.submethod.mode = [];
        tOptions.submethod.toll = [];
        tOptions.submethod.nmax = [];
end


%%% Dimension Check for initial guess
if size(x0,2) > 1
    error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
end


%%% Default method
if isempty(tOptions.method)
    tOptions.method = 'fsolve';      % default implicit solver
end

if isempty(tOptions.theta)
    tOptions.theta = 0.5;            % default theta value
end

                   
%%% Initialization
timerStart = tic;               % timer start
feval = 0;                      % function evaluation counter starts
dimSys = size(x0,1);            % function evaluation step
t = 0:h:tmax;                   % time vector definition
thetaParam = tOptions.theta;
x = [x0 zeros(dimSys,length(t)-1)];    % solution vector allocation


%%% THETA loop
for i = 1 : (length(t)-1)         % main loop of the method
    tk = t(i);
    xk = x(:,i);
    xinter = xk + thetaParam * h * (f(xk,tk));       % x_k+1/2
    
    tn = tk + h;
    thetaFunction = @(xn) xn - xk - h * thetaParam * f(xk,tk) - (1-thetaParam) * h * f(xn,tn) ;
    feval = feval + 2*dimSys;
    
    switch tOptions.method
        case 'fsolve'
            % implicit equation solver
            options = optimoptions ( 'fsolve', 'Display', 'off' );
            [xn,~,conv,info] = fsolve(thetaFunction, xinter, options);    % works in both scalar and vectorial cases

            if conv < 1
                warning('Implicit equation at the step %d has not been solved correctly',i)
            end
            feval = feval + info.funcCount;

        case 'fzero'
            % implicit equation solver
            [xn,~,conv,info] = fzero(thetaFunction, xinter);              % works only in scalar case
            
            if conv < 1
                warning('Implicit equation at the step %d has not been solved correctly',i)
            end
            feval = feval + info.funcCount;

        case 'newton'
            if isempty(tOptions.submethod.mode)
                tOptions.submethod.mode = 'f';        % default implicit solver options
            end
            
            % implicit equation solver
            [xn,conv,info] = newton(thetaFunction, xinter, tOptions.submethod.mode, ...
                                    tOptions.submethod.toll, tOptions.submethod.nmax);    % works only in scalar case for the moment
            if conv == 0 
                warning('Implicit equation at the step %d has not been solved correctly',i)
            end
            feval = feval + info.fevalcost;

        otherwise
            error('Insert a valid method as input\n');
    end  
 
    x(:,i+1) = xn;      % solution value in solution vector
end

elapsedTime = toc(timerStart);   % timer stop

info = struct;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;


%%% plot
if visualConfig.plot == true
    plot(t,x,'o-');
end

end
