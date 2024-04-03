function [x0,infoStartup] = startupGuess(f,x0,order,h,abOptions)
%STARTUPGUESS Summary of this function goes here
%   Detailed explanation goes here
%
%   startup method
%

infoStartup = struct;
    infoStartup.timeCost = 0;
    infoStartup.fevalCost = 0;

visualConf.plot = false;

switch abOptions.startup
    case 'RK'
        % select the runge kutta associated to same order as multistep method
        rkStarter = rkSettings(order,[],[],[]);
        tmaxSingleStep = (order - size(x0,2)) * h;
        [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);
        
    case 'RK1'
        rkStarter = rkSettings(1,[],[],[]);
        tmaxSingleStep = (order - size(x0,2)) * h;
        [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

    case 'RK2'
        rkStarter = rkSettings(2,[],[],[]);
        tmaxSingleStep = (order - size(x0,2)) * h;
        [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);
        
    case 'RK4'
        rkStarter = rkSettings(4,[],[],[]);
        tmaxSingleStep = (order - size(x0,2)) * h;
        [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

    case 'THETA'
        tStarter = tSettings(0.5,[],[]);
        tmaxSingleStep = (order - size(x0,2)) * h;
        [x0,~,infoStartup] = theta(f,x0,tmaxSingleStep,h,tStarter,visualConf);

    case 'IEX4'
        iStarter = iSettings();
        tmaxSingleStep = (order - size(x0,2)) * h;      
        [x0,~,infoStartup] = iex4(f,x0,tmaxSingleStep,h,iStarter,visualConf);

    case 'AB'
        while size(x0,2) ~= order
            abStarter = abSettings(size(x0,2),[],[],[]);
            tmaxLoop = abStarter.order * h;
            [x0,~,infoStartupLoop] = adamsBashforth(f,x0,tmaxLoop,h,abStarter,visualConf);
            
            infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
            infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
        end

    case 'AM'
        while size(x0,2) ~= order
            amStarter = amSettings(size(x0,2),[],[],[],[]);
            tmaxLoop = amStarter.order * h;
            [x0,~,infoStartupLoop] = adamsMoulton(f,x0,tmaxLoop,h,amStarter,visualConf);
            
            infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
            infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
        end
    
    case 'ABM'
        while size(x0,2) ~= order
            abmStarter = abmSettings(size(x0,2),[],[],[],[]);
            tmaxLoop = abmStarter.order * h;
            [x0,~,infoStartupLoop] = adamsBashforthMoulton(f,x0,tmaxLoop,h,abmStarter,visualConf);
            
            infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
            infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
        end
    
    case 'BDF'
        while size(x0,2) ~= order
            bdfStarter = bdfSettings(size(x0,2),[],[],[],[],[]);
            tmaxLoop = bdfStarter.order * h;
            [x0,~,infoStartupLoop] = backwardDifferenceFormula(f,x0,tmaxLoop,h,bdfStarter,visualConf);
            
            infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
            infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
        end

    otherwise
        error('Please insert a valid method as input\n');
end

end