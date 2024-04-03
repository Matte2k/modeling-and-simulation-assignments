function [x,t,info] = adamsMoulton(f,x0,tmax,h,amOptions,visualConfig)

    if nargin < 6
        visualConfig.plot = true;
    end
    
    
    %%% Dimension Check for initial guess
    if size(x0,2) > 1
        warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step\n',h);      % DEBUG
    end
    
    if size(x0,2) > amOptions.order
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
    end
    
    
    %%% Initialization
    infoStartup = [];
    
    
    %%% Adams Bashforth order selector based on 'amOptions'
    switch amOptions.order
        case 1
            [x,t,infoMethod] = am1(f,x0,tmax,h,amOptions);
    
        case 2
            if size(x0,2) < 2
                order = 2;
                [x0,infoStartup] = startupGuess(f,x0,order,h,amOptions);
            end
            
            [x,t,infoMethod] = am2(f,x0,tmax,h,amOptions);
            
        case 3
            if size(x0,2) < 3
                order = 3;
                [x0,infoStartup] = startupGuess(f,x0,order,h,amOptions);
            end
    
            [x,t,infoMethod] = am3(f,x0,tmax,h,amOptions);
    
        case 4
            if size(x0,2) < 4
                order = 4;
                [x0,infoStartup] = startupGuess(f,x0,order,h,amOptions);
            end
    
            [x,t,infoMethod] = am4(f,x0,tmax,h,amOptions);
    
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