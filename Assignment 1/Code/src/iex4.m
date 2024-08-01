function [x,t,info] = iex4(f,x0,tmax,h,iOptions,visualConfig)
%IEX4 - Implicit Extrapolation method application
%
%   Syntax:
%       [x,t,info] = iex4(f,x0,tmax,h,iOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       iOtpions(*),    struct:  see iSettings.m for details
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
%       iOtpions:      set with default solver 'fsolve'
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end

    if nargin < 5
        iOptions.method  = [];
        iOptions.options = [];
    end

    %%% Default submethod
    if isempty(iOptions.method)
        iOptions.method = 'fsolve';      % default implicit solver
    end

    if isempty(iOptions.options)
        iOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );  % default fsolve options
    end

    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected')
    end
                    
    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    alpha = [-1/6, 4, -27/2, 32/3];                 % iex4 parameters definition
    x = [x0 zeros(dimSys,length(t)-1)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];     % fval vector allocation
    feval = feval + dimSys;                         % function evaluation counter update

    %%% IEX4 loop
    for i = 1 : (length(t)-1)         % main loop of the method

        %1st predictor
        tk1 = t(i) + h;
        fk1 = @(k1) x(:,i) + h * f(k1,tk1) - k1;
        [xp1,~,conv(1),info] = fsolve(fk1, x(:,i), iOptions.options);         % initial guess xk
        feval = feval + info.funcCount;

        %2nd predictor
        tk2a = t(i) + h/2;
        fk2a = @(k2a) x(:,i) + h/2 * f(k2a,tk2a) - k2a;
        [k2a,~,conv(2),info] = fsolve(fk2a, x(:,i), iOptions.options);        % initial guess xk
        feval = feval + info.funcCount;

        tk2 = tk1;
        fk2 = @(k2) k2a + h/2 * f(k2,tk2) - k2;
        [xp2,~,conv(3),info] = fsolve(fk2, k2a, iOptions.options);            % initial guess k2a
        feval = feval + info.funcCount;

        %3rd predictor
        tk3a = t(i) + h/3;
        fk3a = @(k3a) x(:,i) + h/3 * f(k3a,tk3a) - k3a;
        [k3a,~,conv(4),info] = fsolve(fk3a, x(:,i), iOptions.options);        % initial guess xk
        feval = feval + info.funcCount;

        tk3b = t(i) + h * 2/3;
        fk3b = @(k3b) k3a + h/3 * f(k3b,tk3b) - k3b;
        [k3b,~,conv(5),info] = fsolve(fk3b, k3a, iOptions.options);           % initial guess k3a
        feval = feval + info.funcCount;

        tk3 = tk1;
        fk3 = @(k3) k3b + h/3 * f(k3,tk3) - k3;
        [xp3,~,conv(6),info] = fsolve(fk3, k3b, iOptions.options);            % initial guess k3b
        feval = feval + info.funcCount;

        %4th predictor
        tk4a = t(i) + h/4;
        fk4a = @(k4a) x(:,i) + h/4 * f(k4a,tk4a) - k4a;
        [k4a,~,conv(7),info] = fsolve(fk4a, x(:,i), iOptions.options);        % initial guess xk
        feval = feval + info.funcCount;

        tk4b = t(i) + h * 2/4;
        fk4b = @(k4b) k4a + h/4 * f(k4b,tk4b) - k4b;
        [k4b,~,conv(8),info] = fsolve(fk4b, k4a, iOptions.options);           % initial guess k4a
        feval = feval + info.funcCount;

        tk4c = t(i) + h * 3/4;
        fk4c = @(k4c) k4b + h/4 * f(k4c,tk4c) - k4c;
        [k4c,~,conv(9),info] = fsolve(fk4c, k4b, iOptions.options);           % initial guess k4b
        feval = feval + info.funcCount;

        tk4 = tk1;
        fk4 = @(k4) k4c + h/4 * f(k4,tk4) - k4;
        [xp4,~,conv(10),info] = fsolve(fk4, k4c, iOptions.options);            % initial guess k4c
        feval = feval + info.funcCount;

        % corrector
        x(:,i+1) = alpha(1) * xp1 + alpha(2) * xp2 + alpha(3) * xp3 + alpha(4) * xp4;
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;     % function evaluation counter update

        % convergence check
        if not(isequal(conv,ones(1,length(conv))))
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end

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