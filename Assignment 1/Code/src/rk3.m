function [x,t,info] = rk3(f,x0,tmax,h,rkOptions)
%RK3 - Runge-Kutta method of order 3
%
%   Syntax:
%       [x,t,info] = rk3(f,x0,tmax,h,rkOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions(*),   struct:  see rkSettings.m for details 
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
%       rkOptions:  set with default alpha and beta
%

    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end

    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'Default';           % Default parameters sets as default
    end

    %%% Parameters definition
    switch rkOptions.method
        case 'Default'
            alpha3 = [1/3 2/3 1]';           
            beta3 = diag([1/3 2/3 3/4 ]);       
            beta3(3,:) = [1/4 0 3/4];

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(rkOptions.alpha) || isempty(rkOptions.beta)            
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(rkOptions.alpha),[3,1])) || not(isequal(size(rkOptions.beta),[3,3]))            
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha3 = rkOptions.alpha;        % custom alpha parameter
            beta3 = rkOptions.beta;          % custom beta parameter

        otherwise
            error('Insert a valid method as input');
    end

    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x =  [x0 zeros(dimSys,length(t)-1)];          % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update

    %%% RK3 loop
    for i=1:(length(t)-1)                              % calculation loop
        fk = fvalVec(:,i);
        xp1 = x(:,i) + beta3(1,1) * h * fk;            % x(i)=xk && t(i)=tk
        tp1 = t(i) + alpha3(1,1) * h;

        xp2 = x(:,i) + beta3(2,2) * h * f(xp1,tp1);    % x(i)=xk && t(i)=tk
        tp2 = t(i) + alpha3(2,1) * h;

        x(:,i+1) = x(:,i) + alpha3(3,1) * h * ( beta3(3,1) * fk ...
            + beta3(3,2) * f(xp1,tp1) ...
            + beta3(3,3) * f(xp2,tp2));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 3*dimSys;             % function evaluation counter update
    end

    elapsedTime = toc(timerStart);   % timer stop

    if nargout == 3
        info = struct;              % info struct build-up
            info.timeCost  = elapsedTime;
            info.fevalCost = feval;
            info.fvalVec   = fvalVec;
            info.implicit  = false;
    end

end