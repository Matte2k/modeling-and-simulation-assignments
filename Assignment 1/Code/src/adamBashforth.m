function [x,t,info] = adamBashforth(f,x0,tmax,h,abOptions,visualConfig)

if nargin < 6
    visualConfig.plot = true;
end

%%% Runge-Kutta method selector based on 'rkOptions'
switch rkOptions.method
    case 'AB1'
        [x,t,info] = ab1(f,x0,tmax,h,abOptions);

    case 'AB2'
        [x,t,info] = ab2(f,x0,tmax,h,abOptions);

    case 'AB3'
        [x,t,info] = ab4(f,x0,tmax,h,abOptions);

    case 'AB4'
        [x,t,info] = ab4(f,x0,tmax,h,abOptions);

    otherwise
        error('Please insert a valid method as input\n');
end

if visualConfig.plot == true
    plot(t,x,'o-');
end

end