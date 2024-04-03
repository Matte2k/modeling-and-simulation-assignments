function [x,t,info] = backwardDifferenceFormula(f,x0,tmax,h,bdfOptions,visualConfig)

    if nargin < 6
        visualConfig.plot = true;
    end
    
    
    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step\n',h);      % DEBUG
    end
    
    if size(x0,2) > bdfOptions.order
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
    end
    
    
    %%% Initialization
    infoStartup = [];
    
    
    %%% Adams Bashforth Moulton order selector based on 'bdfOptions'
    switch bdfOptions.order
        case 1
            [x,t,infoMethod] = bdf1(f,x0,tmax,h,bdfOptions);
    
        case 2
            if size(x0,2) < 2
                order = 2;
                [x0,infoStartup] = startupGuess(f,x0,order,h,bdfOptions);
            end
            
            [x,t,infoMethod] = bdf2(f,x0,tmax,h,bdfOptions);
            
        case 3
            if size(x0,2) < 3
                order = 3;
                [x0,infoStartup] = startupGuess(f,x0,order,h,bdfOptions);
            end
    
            [x,t,infoMethod] = bdf3(f,x0,tmax,h,bdfOptions);
    
        case 4
            if size(x0,2) < 4
                order = 4;
                [x0,infoStartup] = startupGuess(f,x0,order,h,bdfOptions);
            end
    
            [x,t,infoMethod] = bdf4(f,x0,tmax,h,bdfOptions);
    
        otherwise
            error('Please insert a valid order as input\n');
    end
    
    if visualConfig.plot == true
        plot(t,x,'o-');
    end
    
    
    if nargout == 3
        info = struct;
    
        % Merge method and startup info (if present)
        if not(isempty(infoStartup))
            info.timeCost = infoMethod.timeCost + infoStartup.timeCost;
            info.fevalCost = infoMethod.fevalCost + infoStartup.fevalCost;
        else
            info.timeCost = infoMethod.timeCost;
            info.fevalCost = infoMethod.fevalCost;
        end
    end
    
    
    end