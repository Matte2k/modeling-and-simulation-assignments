function [x,t,info] = bdf1(f,x0,tmax,h,bdfOptions)
%BDF1 - Backward Difference Formula method of order 1
%
%   Syntax:
%       [x,t,info] = bdf1(f,x0,tmax,h,bdfOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       bdfOtpions(*),  struct:  see bdfSettings.m for details 
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
%       bdfOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        bdfOptions.method  = [];
        bdfOptions.options = [];
        bdfOptions.alpha   = [];
        bdfOptions.beta    = [];
        bdfOptions.gamma   = [];
    end


    %%% Default method
    if isempty(bdfOptions.method)
        bdfOptions.method = 'Standard';
    end

    if isempty(bdfOptions.options)
        bdfOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );   % default fsolve options
    end


    %%% Parameters definition
    switch bdfOptions.method
        case 'Standard'
            alpha = 1;      % standard alpha parameter
            beta = 1;       % standard beta parameter
            gamma = 1;      % standard gamma parameter

            if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))            % TO DEBUG
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    bdfOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)            % TO DEBUG
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,1]))            % TO DEBUG
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = bdfOptions.alpha;   % custom alpha parameter
            beta = bdfOptions.beta;     % custom beta parameter
            gamma = bdfOptions.gamma;   % custom gamma parameter

        otherwise
            error('Insert a valid method as input');

    end


    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)]; % fval vector allocation
    feval = feval + dimSys;                     % function evaluation counter update


    %%% BDF1 loop
    for i = 1 : (length(t)-1)       % main loop of the method
        xk1 = x(:,i);
        tp = t(i) + h;
        fp = @(xp) 1/alpha * ( gamma(1)*xk1 + beta*h*f(xp,tp) ) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), bdfOptions.options);
        feval = feval + info.funcCount;

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

end