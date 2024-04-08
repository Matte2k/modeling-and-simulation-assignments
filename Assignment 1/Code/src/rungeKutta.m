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
            [x,t,info] = rk1(f,x0,tmax,h,rkOptions);

        case 2
            [x,t,info] = rk2(f,x0,tmax,h,rkOptions);

        case 3
            [x,t,info] = rk3(f,x0,tmax,h,rkOptions);

        case 4
            [x,t,info] = rk4(f,x0,tmax,h,rkOptions);

        otherwise
            error('Please insert a valid method as input\n');
    end

    if visualConfig == true
        plot(t,x,'o-');
    end

end