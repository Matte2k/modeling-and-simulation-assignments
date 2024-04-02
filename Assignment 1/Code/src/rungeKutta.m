function [x,t,info] = rungeKutta(f,x0,tmax,h,rkOptions,visualConfig)
%RUNGEKUTTA Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    visualConfig.plot = true;
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
        error('TO BE ADDED\n');

    case 4
        [x,t,info] = rk4(f,x0,tmax,h,rkOptions);

    otherwise
        error('Please insert a valid method as input\n');
end

if visualConfig.plot == true
    plot(t,x,'o-');
end

end

