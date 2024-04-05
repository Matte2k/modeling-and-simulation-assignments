clearvars; close all; clc

f = @(x,t) [-5/2.*(1+8*sin(t)).*x(1); (1-x(1)).*x(2) + x(1)];
x0 = [1 1]';
tmax = 3;
h = 0.1;

options = optimoptions ( 'fsolve', 'Display', 'off' );

m = length ( x0 );              % sys dimension
t = 0:h:tmax;
x = zeros ( m , length(t) );    % solution vector
x(:,1) = x0(:);                 % initial value definition

%%% BDF
%   implement generic BDF method in order to have the possiblity to use 2
%   different startup

alpha = [ 1   3   11  25];

beta =  [ 1   2   6   12];

gamma =  [ 1     0    0   0; ...
           4    -1    0   0; ...
           18   -9    2   0; ...
           48   -36   16  3];

%BDF1
for i = 1 : (length(t)-1)         % main loop of the method
    xk1 = x(:,i);
    
    tp = t(i) + h;
    fp = @(xp) 1/alpha(1) * ( gamma(1,1)*xk1 + beta(1)*h*f(xp,tp) ) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk

end


%BDF2
for i = 2 : (length(t)-1)         % main loop of the method
    xk1 = x(:,i  );
    xk2 = x(:,i-1);
    
    tp = t(i) + h;
    fp = @(xp) 1/alpha(2) * ( gamma(2,1)*xk1 + gamma(2,2)*xk2 + beta(2)*h*f(xp,tp) ) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end


%BDF3
for i = 3 : (length(t)-1)         % main loop of the method
    xk1 = x(:,i  );
    xk2 = x(:,i-1);
    xk3 = x(:,i-2);

    tp = t(i) + h;
    fp = @(xp) 1/alpha(3) * ( gamma(3,1)*xk1 + gamma(3,2)*xk2 + gamma(3,3)*xk3 + beta(3)*h*f(xp,tp) ) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end


%BDF4
for i = 4 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    fk4 = f(x(:,i-2),t(i-3));

    tp = t(i) + h;
    fp = @(xp) x(:,i) + h/alpha(4) * (beta(4,1)*f(xp,tp) + beta(4,2)*fk1 + beta(4,3)*fk2 + beta(4,4)*fk3) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end


% TO ADD: startup check for BDFn