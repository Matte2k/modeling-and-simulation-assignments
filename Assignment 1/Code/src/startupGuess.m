function [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,methodOptions)
%STARTUP GUESS - Solve the startup problem for multistep method
%
%   Syntax:
%       [x0,infoStartup] = startupGuess(f,x0,guessOrder,h,methodOptions)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       guessOrder,     double:  number of element in initial guess
%       h,              double:  time step of the integration
%       methodOptions,  struct:  see corresponding method Settings.m for details
%
%   Output:
%       x,     double[n,m]:  solution vector
%       info,       struct:  information on method used:
%           - info.timeCost,     double:  time spent
%           - info.fevalCost,    double:  # of function evaluations
%           - info.fvalVec, dobule[n,m]:  f evaluated in solution points
%           - info.implicit,       bool:  true if the method is implicit
%


    %%% Initialization
    if nargout == 2
        infoStartup = struct;
            infoStartup.timeCost  = [];
            infoStartup.fevalCost = [];
            infoStartup.implicit  = [];
    end
    visualConf = false;     % disable plot in the startup problem


    %%% Startup solver selector vased on 'methodOptions'
    switch methodOptions.startup
        case 'RK'
            % select the runge kutta associated to same order as multistep method
            rkStarter = rkSettings(guessOrder);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);
            
        case 'RK1'
            rkStarter = rkSettings(1);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

        case 'RK2'
            rkStarter = rkSettings(2);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

        case 'RK3'
            rkStarter = rkSettings(3);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);
            
        case 'RK4'
            rkStarter = rkSettings(4);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = rungeKutta(f,x0,tmaxSingleStep,h,rkStarter,visualConf);

        case 'THETA'
            tStarter = tSettings(0.5);
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;
            [x0,~,infoStartup] = theta(f,x0,tmaxSingleStep,h,tStarter,visualConf);

        case 'IEX4'
            iStarter = iSettings();
            tmaxSingleStep = (guessOrder - size(x0,2)) * h;      
            [x0,~,infoStartup] = iex4(f,x0,tmaxSingleStep,h,iStarter,visualConf);

        case 'AB'
            while size(x0,2) ~= guessOrder
                abStarter = abSettings(size(x0,2));
                tmaxLoop = abStarter.order * h;
                [x0,~,infoStartupLoop] = adamsBashforth(f,x0,tmaxLoop,h,abStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end

        case 'AM'
            while size(x0,2) ~= guessOrder
                amStarter = amSettings(size(x0,2));
                tmaxLoop = amStarter.order * h;
                [x0,~,infoStartupLoop] = adamsMoulton(f,x0,tmaxLoop,h,amStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end
        
        case 'ABM'
            while size(x0,2) ~= guessOrder
                abmStarter = abmSettings(size(x0,2));
                tmaxLoop = abmStarter.order * h;
                [x0,~,infoStartupLoop] = adamsBashforthMoulton(f,x0,tmaxLoop,h,abmStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end
        
        case 'BDF'
            while size(x0,2) ~= guessOrder
                bdfStarter = bdfSettings(size(x0,2));
                tmaxLoop = bdfStarter.order * h;
                [x0,~,infoStartupLoop] = backwardDifferenceFormula(f,x0,tmaxLoop,h,bdfStarter,visualConf);
                
                infoStartup.timeCost = infoStartup.timeCost + infoStartupLoop.timeCost;
                infoStartup.fevalCost = infoStartup.fevalCost + infoStartupLoop.fevalCost;
            end

        otherwise
            error('Please insert a valid method as input');
    end

end