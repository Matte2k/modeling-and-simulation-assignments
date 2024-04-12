function [x,t,info] = theta(f,x0,tmax,h,tOptions,visualConfig)
%THETA - Theta method application
%
%   Syntax:
%       [x,t,info] = theta(f,x0,tmax,h,tOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       tOtpions(*),    struct:  see tSettings.m for details
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       x,     double[n,m]:  solution vector
%       t,     double[1,m]:  time istant associated to solutions
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%
%   Default settings for optional input (*):
%       tOtpions:      set with default solver 'fsolve' and 'theta = 0.5'
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    if nargin < 5
    tOptions = struct;
        tOptions.theta  = [];
        tOptions.method = [];
    end


    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected')
    end

    if size(x0,1) > 1 && not(isequal(tOptions.method,'fsolve'))
        error('The method select currently do not support multivariable system')
    end


    %%% Default method
    if isempty(tOptions.theta)
        tOptions.theta = 0.5;            % default theta value
    end

    if isempty(tOptions.method)
        tOptions.method = 'fsolve';      % default implicit solver
    end

                    
    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    thetaParam = tOptions.theta;
    x = [x0 zeros(dimSys,length(t)-1)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];     % fval vector allocation
    feval = feval + dimSys;                         % function evaluation counter update


    %%% THETA loop
    for i = 1 : (length(t)-1)                   % main loop of the method
        tk = t(i);
        xk = x(:,i);
        fk = fvalVec(:,i);                      % fk = f(xk,tk)
        xinter = xk + thetaParam * h * fk;      % x_k+1/2
        
        tn = tk + h;
        thetaFunction = @(xn) xn - xk - h * thetaParam * fk - (1-thetaParam) * h * f(xn,tn) ;
        
        switch tOptions.method
            case 'fsolve'
                % implicit equation solver
                options = optimoptions ( 'fsolve', 'Display', 'off' );
                [xn,fvalVec(:,i+1),conv,info] = fsolve(thetaFunction, xinter, options);      % works in both scalar and vectorial cases

                if conv < 1
                    warning('Implicit equation at the step %d has not been solved correctly',i)
                end
                feval = feval + info.funcCount;

            case 'fzero'
                % implicit equation solver
                [xn,fvalVec(:,i+1),conv,info] = fzero(thetaFunction, xinter);                % works only in scalar case
                
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
                                        tOptions.submethod.toll, tOptions.submethod.nmax);   % works only in scalar case for the moment
                if conv == 0 
                    warning('Implicit equation at the step %d has not been solved correctly',i)
                end
                fvalVec(:,i+1) = f(xn,t(i+1));
                feval = feval + info.fevalcost + dimSys;

            otherwise
                error('Insert a valid method as input');
        end  
    
        x(:,i+1) = xn;              % solution value in solution vector
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = true;
    end


    %%% plot
    if visualConfig == true
        plot(t,x,'o-');
    end

end