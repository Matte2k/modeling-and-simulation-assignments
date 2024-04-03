function [x,t,info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptions,visualConfig)

if nargin < 6
    visualConfig.plot = true;
end


%%% Dimension Check for initial guess
if size(x0,2) > 1
    warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step\n',h);      % DEBUG
end

if size(x0,2) > abmOptions.order
    error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
end


%%% Initialization
infoStartup = [];


%%% Adams Bashforth Moulton order selector based on 'abmOptions'
switch abmOptions.order
    case 1
        [x,t,infoMethod] = abm1(f,x0,tmax,h,abmOptions);

    case 2
        if size(x0,2) < 2
            order = 2;
            [x0,infoStartup] = startupGuess(f,x0,order,h,abmOptions);
        end
        
        [x,t,infoMethod] = abm2(f,x0,tmax,h,abmOptions);
        
    case 3
        if size(x0,2) < 3
            order = 3;
            [x0,infoStartup] = startupGuess(f,x0,order,h,abmOptions);
        end

        [x,t,infoMethod] = abm3(f,x0,tmax,h,abmOptions);

    case 4
        if size(x0,2) < 4
            order = 4;
            [x0,infoStartup] = startupGuess(f,x0,order,h,abmOptions);
        end

        [x,t,infoMethod] = abm4(f,x0,tmax,h,abmOptions);

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