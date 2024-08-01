function [x,t,info] = ab3(f,x0,tmax,h,abOptions)
%AB3 - Adams Bashforth method of order 3
%
%   Syntax:
%       [x,t,info] = ab3(f,x0,tmax,h,abOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,3]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abOtpions(*),   struct:  see abSettings.m for details 
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
%       abOptions:  set with default alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        abOptions.method = [];
        abOptions.alpha  = [];
        abOptions.beta   = [];
    end

    %%% Default method
    if isempty(abOptions.method)
        abOptions.method = 'Standard';
    end

    %%% Parameters definition
    switch abOptions.method
        case 'Standard'
            alpha = 12;             % standard alpha parameter
            beta = [23 -16  5];     % standard beta parameter

            if not(isempty(abOptions.alpha)) || not(isempty(abOptions.beta))
                warning('Parameters matrix unused, %s parameters are used instead',...
                    abOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(abOptions.alpha) || isempty(abOptions.beta)
                error('No parameters matrix has been given as input');    % missing parameters matrix

            elseif not(isequal(size(abOptions.alpha),[1,1])) || not(isequal(size(abOptions.beta),[1,3]))
                error('Parameters matrix dimensions are invalid');        % parameters matrix with wrong size

            end

            alpha = abOptions.alpha;    % custom alpha parameter
            beta = abOptions.beta;      % custom beta parameter

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

    %%% AB3 loop
    for i = 3 : (length(t)-1)        % main loop of the method
        fk1 = fvalVec(:,i);
        fk2 = fvalVec(:,i-1);
        fk3 = fvalVec(:,i-2);
        x(:,i+1) = x(:,i) + h/alpha * (beta(1)*fk1 + beta(2)*fk2 + beta(3)*fk3);
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;      % function evaluation counter updateend
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;               % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end