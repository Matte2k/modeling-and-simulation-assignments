function [x,t,info] = abm3(f,x0,tmax,h,abmOptions)
%ABM3 - Adams Bashforth Moulton method of order 3
%
%   Syntax:
%       [x,t,info] = abm3(f,x0,tmax,h,abmOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,3]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions(*),  struct:  see abmSettings.m for details 
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
%       abmOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abmOptions.method = [];
        abmOptions.alpha  = [];
        abmOptions.betaP  = [];
        abmOptions.betaC  = [];
    end

    %%% Default method
    if isempty(abmOptions.method)
        abmOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abmOptions.method
        case 'Standard'
            alpha = 12;             % standard alpha parameter
            betaP = [23 -16  5];    % standard beta predictor parameter
            betaC = [5   8  -1];    % standard beta corrector parameter

            if not(isempty(abmOptions.alpha)) || not(isempty(abmOptions.betaP)) || not(isempty(abmOptions.betaC))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    abmOptions.method);     % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abmOptions.alpha) || isempty(abmOptions.betaP) || isempty(abmOptions.betaC)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(abmOptions.alpha),[1,1])) || not(isequal(size(abmOptions.betaP),[1,3])) || not(isequal(size(abmOptions.betaC),[1,3]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = abmOptions.alpha;     % custom alpha parameter
            betaP = abmOptions.betaP;     % custom beta predictor parameter
            betaC = abmOptions.betaC;     % beta corrector parameter

        otherwise
            error('Insert a valid method as input');

    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-3)];             % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), f(x(:,3),t(3)), ...
                    zeros(dimSys,length(t)-3)];     % fval vector allocation
    feval = feval + 3*dimSys;                       % function evaluation counter update

    %%% ABM3 loop
    for i = 3 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        xp = x(:,i) + h/alpha * (betaP(1)*fk1 + betaP(2)*fk2 + betaP(3)*fk3);
        tp = t(i) + h;
        x(:,i+1) = x(:,i) + h/alpha * (betaC(1)*f(xp,tp) + betaC(2)*fk1 + betaC(3)*fk2);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;   % function evaluation counter updateend
    end

    elapsedTime = toc(timerStart);  % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end