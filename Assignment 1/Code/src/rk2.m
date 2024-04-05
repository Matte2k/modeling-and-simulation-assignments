function [x,t,info] = rk2(f,x0,tmax,h,rkOptions)
%RK2 - Runge-Kutta method of order 2
%
%   Syntax:
%       [x,t,info] = rk2(f,x0,tmax,h,rkOptions)
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
%       rkOptions:  set with default 'Heun' alpha and beta
%


    %%% Optional input definition
    if nargin < 5
        rkOptions.method = [];
        rkOptions.alpha  = [];
        rkOptions.beta   = [];
    end


    %%% Default method
    if isempty(rkOptions.method)
        rkOptions.method = 'Heun';          % Heun parameters sets as default
    end


    %%% Parameters definition
    switch rkOptions.method
        case 'Heun'
            alpha2 = [1 1]';                 % rk alpha matrix
            beta2 = [1, 0; 0.5, 0.5];        % rk beta matrix

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'MidPoint'
            alpha2 = [0.5 1]';               % rk alpha matrix
            beta2 = [0.5, 0; 0, 1];          % rk beta matrix

            if not(isempty(rkOptions.alpha)) || not(isempty(rkOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    rkOptions.method);       % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(rkOptions.alpha) || isempty(rkOptions.beta)
                error('No parameters matrix has been given as input')   % missing parameters matrix

            elseif not(isequal(size(rkOptions.alpha),[2,1])) || not(isequal(size(rkOptions.beta),[2,2]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size
                
            end

            alpha2 = rkOptions.alpha;        % custom alpha parameter
            beta2 = rkOptions.beta;          % custom beta parameter

            %%% Check if the custom parameters are valid
            eC1 = beta2(2,1) + beta2(2,2) == 1;
            eC2 = 2 * beta2(1,1) * beta2(2,2) == 1;
            eC3 = 2 * alpha2(1,1) * beta2(2,2) == 1;

            if eC1 ~= true ||  eC2 ~= true ||  eC3 ~= true
                error('Insert a valid parameters matrix as input');
            end

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


    %%% RK2 loop
    for i=1:(length(t)-1)
        fk = fvalVec(:,i);
        xp = x(:,i) + beta2(1,1) * h * fk;    % x(i)=xk && t(i)=tk
        tp = t(i) + alpha2(1,1) * h;
        x(:,i+1) = x(:,i) + alpha2(2,1) * h * (beta2(2,1) * fk + beta2(2,2) * f(xp,tp));
        fvalVec(:,i+1) = f(x(:,i+1),t(i+1));
        feval = feval + 2*dimSys;             % function evaluation counter update
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