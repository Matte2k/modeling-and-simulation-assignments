function [x,t,info] = adamsBashforth(f,x0,tmax,h,abOptions,visualConfig)

if nargin < 6
    visualConfig.plot = true;
end


%%% Dimension Check for initial guess
if size(x0,2) > 1
    warning('INPUT initial guess MUST be evaluated on an equally spaced grid with %.3d as step\n',h);      % DEBUG
end

if size(x0,2) > abOptions.order
    error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
end


%%% Initialization
infoStartup = [];


%%% Adams Bashforth order selector based on 'abOptions'
switch abOptions.order
    case 1
        [x,t,infoMethod] = ab1(f,x0,tmax,h,abOptions);

    case 2
        if size(x0,2) < 2
            order = 2;
            [x0,infoStartup] = startupGuess(f,x0,order,h,abOptions);
        end
        
        [x,t,infoMethod] = ab2(f,x0,tmax,h,abOptions);
        
    case 3
        if size(x0,2) < 3
            order = 3;
            [x0,infoStartup] = startupGuess(f,x0,order,h,abOptions);
        end

        [x,t,infoMethod] = ab3(f,x0,tmax,h,abOptions);

    case 4
        if size(x0,2) < 4
            order = 4;
            [x0,infoStartup] = startupGuess(f,x0,order,h,abOptions);
        end

        [x,t,infoMethod] = ab4(f,x0,tmax,h,abOptions);

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