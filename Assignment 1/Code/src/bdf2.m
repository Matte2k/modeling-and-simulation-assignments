function [x,t,info] = bdf2(f,x0,tmax,h,bdfOptions)
%BDF2 - Backward Difference Formula method of order 2
%
%   Syntax:
%       [x,t,info] = bdf2(f,x0,tmax,h,bdfOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,2]:  inital guess 
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
            alpha = 3;          % standard alpha parameter
            beta = 2;           % standard beta parameter
            gamma = [4 -1];     % standard gamma parameter
            
            if not(isempty(bdfOptions.alpha)) || not(isempty(bdfOptions.beta)) || not(isempty(bdfOptions.gamma))
                warning('Parameters matrix unused, standard %s parameters are used instead',...
                    bdfOptions.method);     % warning for parameters matrix unused
            end
            
        case 'Custom'
            if isempty(bdfOptions.alpha) || isempty(bdfOptions.beta) || isempty(bdfOptions.gamma)
                error('No parameters matrix has been given as input\n');  % missing parameters matrix
                
            elseif not(isequal(size(bdfOptions.alpha),[1,1])) || not(isequal(size(bdfOptions.beta),[1,1])) || not(isequal(size(bdfOptions.gamma),[1,2]))
                error('Parameters matrix dimensions are invalid\n');      % parameters matrix with wrong size
                
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

    x = [x0 zeros(dimSys,length(t)-2)];         % solution vector allocation
    fvalVec = [f(x(:,1),t(1)), f(x(:,2),t(2)), ...
                    zeros(dimSys,length(t)-2)]; % fval vector allocation
    feval = feval + 2*dimSys;                   % function evaluation counter update


    %%% BDF2 loop
    for i = 2 : (length(t)-1)       % main loop of the method
        xk1 = x(:,i  );
        xk2 = x(:,i-1);
        tp = t(i) + h;
        fp = @(xp) 1/alpha * ( gamma(1)*xk1 + gamma(2)*xk2 + beta*h*f(xp,tp) ) - xp;
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