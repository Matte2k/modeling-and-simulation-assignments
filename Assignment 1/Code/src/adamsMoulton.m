function [x,t,info] = adamsMoulton(f,x0,tmax,h,amOptions,visualConfig)
%ADAMS MOULTON - Adams Moulton method selection and application
%
%   Syntax:
%       [x,t,info] = adamsMoulton(f,x0,tmax,h,amOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       amOtpions,      struct:  see amSettings.m for details
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
%       visualConfig:  set as true by default
%


    %%% Optional input definition
    if nargin < 6
        visualConfig = true;
    end
        
        
    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        % warning to remember that initial guess correct input format
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step',h);      % DEBUG
    end

    % case for AM order grater than 1
    if amOptions.order > 1
        if size(x0,2) > amOptions.order-1 
            error('The initial guess is invalid, too many input in x0 vector compare to the order selected')

        elseif size(x0,2) < amOptions.order-1 
            if isempty(amOptions.startup)
                error('Need to define startup method to retrive a valid initial guess')
            end

        else
            infoStartup = [];   % initialization of the infoStartup variable used in info
            if not(isempty(amOptions.startup))
                warning('Startup method defined will not be used, x0 is already valid')
            end
        end
    end

    % case for AM order equal to 1
    if amOptions.order == 1
        if size(x0,2) > amOptions.order 
            error('The initial guess is invalid, too many input in x0 vector compare to the order selected')
        
        elseif size(x0,2) < amOptions.order 
            if isempty(amOptions.startup)
                error('Need to define startup method to retrive a valid initial guess')
            end
        
        else
            infoStartup = [];   % initialization of the infoStartup variable used in info
            if not(isempty(amOptions.startup))
                warning('Startup method defined will not be used, x0 is already valid')
            end
        end
    end


    %%% Adams Bashforth order selector based on 'amOptions'
    switch amOptions.order
        case 1
            [x,t,infoMethod] = am1(f,x0,tmax,h,amOptions);

        case 2
            [x,t,infoMethod] = am2(f,x0,tmax,h,amOptions);
            
        case 3
            % retrive a valid initial guess based on 'amOptions'
            if size(x0,2) < 2
                guessOrder = 2;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,amOptions);
            end

            [x,t,infoMethod] = am3(f,x0,tmax,h,amOptions);

        case 4
            % retrive a valid initial guess based on 'amOptions'
            if size(x0,2) < 3
                guessOrder = 3;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,amOptions);
            end

            [x,t,infoMethod] = am4(f,x0,tmax,h,amOptions);

        otherwise
            error('Please insert a valid order as input');
    end

    if visualConfig == true
        plot(t,x,'o-');     % plot of the solution
    end


    if nargout == 3
        info = struct;
            info.fvalVec  = infoMethod.fvalVec; 
            info.implicit = true;

        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost  = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
        else
            info.timeCost  = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
        end
    end

end