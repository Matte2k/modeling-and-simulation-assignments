function [x,t,info] = rk1(f,x0,tmax,h,rkOptions)
%RK1 - Runge-Kutta method of order 1
%
%   Syntax:
%       [x,t,info] = rk1(f,x0,tmax,h,rkOptions)
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
%       rkOptions:  set with default 'FowardEuler' alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end


    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'FowardEuler';   % FowardEuler parameters sets as default
    end

    if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))            % TO DEBUG
        warning('Parameters matrix unused, standard %s parameters are used instead',...
            rkOptions.method);      % warning for parameters matrix unused
    end


    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = size(x0,1);            % function evaluation step
    t = 0:h:tmax;                   % time vector definition

    x = [x0 zeros(dimSys,length(t)-1)];           % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), ...
                    zeros(dimSys,length(t)-1)];   % fval vector allocation
    feval = feval + dimSys;                       % function evaluation counter update


    %%% RK1 loop
    for i=1:(length(t)-1)
        x(:,i+1) = x(:,i) + h * f(x(:,i),t(i));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + dimSys;      % function evaluation counter update
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