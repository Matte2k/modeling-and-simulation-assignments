function [x,t,info] = rungeKutta(f,x0,tmax,h,rkOptions,visualConfig)
%RUNGEKUTTA Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    visualConfig.plot = true;
end

%%% Runge-Kutta method selector based on 'rkOptions'
switch rkOptions.method
    case 'RK1'
        [x,t,info] = rk1(f,x0,tmax,h,rkOptions);

    case 'RK2'
        [x,t,info] = rk2(f,x0,tmax,h,rkOptions);

    case 'RK4'
        [x,t,info] = rk4(f,x0,tmax,h,rkOptions);

    otherwise
        error('Please insert a valid method as input\n');
end

if visualConfig.plot == true
    plot(t,x,'o-');
end

end

