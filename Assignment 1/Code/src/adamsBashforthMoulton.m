function [x,t,info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptions,visualConfig)
%ADAMS BASHFORTH MOULTON - Adams Bashforth Moulton method selection and application
%
%   Syntax:
%       [x,t,info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       abmOtpions,      struct:  see abmSettings.m for details
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

    if size(x0,2) > abmOptions.order
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')

    elseif size(x0,2) < abmOptions.order 
        if isempty(abmOptions.startup)
            error('Need to define startup method to retrive a valid initial guess')
        end

    else
        infoStartup = [];   % initialization of the infoStartup variable used in info
        if not(isempty(abmOptions.startup))
            warning('Startup method defined will not be used, x0 is already valid')
        end
    end


    %%% Adams Bashforth Moulton order selector based on 'abmOptions'
    switch abmOptions.order
        case 1
            [x,t,infoMethod] = abm1(f,x0,tmax,h,abmOptions);

        case 2
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 2
                guessOrder = 2;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abmOptions);
            end
            
            [x,t,infoMethod] = abm2(f,x0,tmax,h,abmOptions);
            
        case 3
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 3
                guessOrder = 3;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abmOptions);
            end

            [x,t,infoMethod] = abm3(f,x0,tmax,h,abmOptions);

        case 4
            % retrive a valid initial guess based on 'abOptions'
            if size(x0,2) < 4
                guessOrder = 4;
                [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,abmOptions);
            end

            [x,t,infoMethod] = abm4(f,x0,tmax,h,abmOptions);

        otherwise
            error('Please insert a valid order as input');
    end

    if visualConfig == true
        plot(t,x,'o-');     % plot of the solution
    end


    if nargout == 3
        info = struct;      % info struct build-up
            info.fvalVec   = infoMethod.fvalVec;        

        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost  = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
            
            % update the implict flag checking the startup method implicity
            if infoStartup.implicit == true
                info.implicit = true;
            else
                info.implicit = false;
            end
        else
            info.timeCost  = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
            info.implicit  = false;
        end
    end

end