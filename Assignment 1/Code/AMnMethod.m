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

% TO DO: ABM3, BDF3

%%% ABM3
%   implement generic AM method in order to have the possiblity to use 2
%   different startup

alpha = [ 1   2   12  24];

beta =  [ 1   0   0   0; ...
          1   1   0   0; ...
          5   8   -1  0; ...
          9   19  -5  1];

%AM1
for i = 1 : (length(t)-1)         % main loop of the method
    tp = t(i) + h;
    fp = @(xp) x(:,i) + h/alpha(1) * (beta(1,1)*f(xp,tp)) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end

%AM2 - crank nicholson
for i = 2 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    
    tp = t(i) + h;
    fp = @(xp) x(:,i) + h/alpha(2) * (beta(2,1)*f(xp,tp) + beta(2,2)*fk1) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end

%AM3
for i = 3 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    
    tp = t(i) + h;
    fp = @(xp) x(:,i) + h/alpha(3) * (beta(3,1)*f(xp,tp) + beta(3,2)*fk1 + beta(3,3)*fk2) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end

%AM4
for i = 4 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    
    tp = t(i) + h;
    fp = @(xp) x(:,i) + h/alpha(4) * (beta(4,1)*f(xp,tp) + beta(4,2)*fk1 + beta(4,3)*fk2 + beta(4,4)*fk3) - xp;   % k+1
    x(:,i+1) = fsolve(fp, x(:,i), options);         % initial guess xk
end

% TO ADD: startup check for AMn 