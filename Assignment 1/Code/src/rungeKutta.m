function [x,t,info] = rungeKutta(f,x0,tmax,h,rkOptions,visualConfig)
%RUNGE KUTTA - Runge-Kutta method selection and application
%
%   Syntax:
%       [x,t,info] = rungeKutta(f,x0,tmax,h,rkOptions,visualConfig)
%
%   Input:
%       f,       function(x,t):  IVP problem
%       x0,        double[n,#]:  generic initial guess 
%       tmax,           double:  upper time limit of the integration
%       h,              double:  time step of the integration
%       rkOtpions,      struct:  see rkSettings.m for details
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
%           - info.iterations,    double:
%           - info.avgTimeCost,  double:  averge time spent on 'avgTime' iterations
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
        error('The initial guess is invalid, too many input in x0 vector compare to the order selected\n')
    end

    
    %%% Runge-Kutta method selector based on 'rkOptions'
    switch rkOptions.order
        case 1 
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk1(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        case 2
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk2(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        case 3
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk3(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        case 4
            tic
            for i = 1:rkOptions.iterations
                [x,t,info] = rk4(f,x0,tmax,h,rkOptions);
            end
            info.avgTimeCost = toc/rkOptions.iterations;

        otherwise
            error('Please insert a valid method as input\n');
    end
    info.iterations = rkOptions.iterations;


    if visualConfig == true
        plot(t,x,'-');
    end

end