function [x,t,info] = am2(f,x0,tmax,h,amOptions)
%AM2 - Adams Moulton method of order 2
%
%   Syntax:
%       [x,t,info] = am2(f,x0,tmax,h,amOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,1]:  inital guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions(*),   struct:  see amSettings.m for details 
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
%       amOptions:  set with default alpha, beta and solver options
%


    %%% Optional input definition
    if nargin < 5
        amOptions.method  = [];
        amOptions.options = [];
        amOptions.alpha   = [];
        amOptions.beta    = [];
    end


    %%% Default method
    if isempty(amOptions.method)
        amOptions.method = 'Standard';
    end

    if isempty(amOptions.options)
        amOptions.options = optimoptions ( 'fsolve', 'Display', 'off' );    % default fsolve options
    end


    %%% Parameters definition
    switch amOptions.method
        case 'Standard'
            alpha = 2;      % standard alpha parameter
            beta = [1 1];   % standard beta parameter

            if not(isempty(amOptions.alpha)) || not(isempty(amOptions.beta))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    amOptions.method);      % warning for parameters matrix unused
            end

        case 'Custom'
            if isempty(amOptions.alpha) || isempty(amOptions.beta)
                error('No parameters matrix has been given as input');  % missing parameters matrix

            elseif not(isequal(size(amOptions.alpha),[1,1])) || not(isequal(size(amOptions.beta),[1,2]))
                error('Parameters matrix dimensions are invalid');      % parameters matrix with wrong size

            end

            alpha = amOptions.alpha;    % custom alpha parameter
            beta = amOptions.beta;      % custom beta parameter

        otherwise
            error('Insert a valid method as input');

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


    %%% AM2 loop
    for i = 1 : (length(t)-1)       % main loop of the method
        fk1 = fvalVec(:,i);
        tp = t(i) + h;
        fp = @(xp) x(:,i) + h/alpha * (beta(1)*f(xp,tp) + beta(2)*fk1) - xp;
        [x(:,i+1),fvalVec(:,i+1),conv,info] = fsolve(fp, x(:,i), amOptions.options);
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