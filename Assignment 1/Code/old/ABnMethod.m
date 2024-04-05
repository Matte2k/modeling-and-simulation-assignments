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

%%% AB3
%   implement generic AB method in order to have the possiblity to use 2
%   different startup

alpha = [ 1   2   12  24];

beta =  [ 1   0   0   0; ...
          3  -1   0   0; ...
          23 -16  5   0; ...
          55 -59  37 -9];

%AB1
for i = 1 : (length(t)-1)         % main loop of the method
    fk = f(x(:,i),t(i));
    x(:,i+1) = x(:,i) + h/alpha(1) * (beta(1)*fk);
end

%AB2
for i = 2 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    x(:,i+1) = x(:,i) + h/alpha(2) * (beta(2,1)*fk1 + beta(2,2)*fk2);
end

%AB3
for i = 3 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    x(:,i+1) = x(:,i) + h/alpha(3) * (beta(3,1)*fk1 + beta(3,2)*fk2 + beta(3,3)*fk3);
end

%AB4
for i = 4 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    fk4 = f(x(:,i-3),t(i-3));
    x(:,i+1) = x(:,i) + h/alpha(4) * (beta(4,1)*fk1 + beta(4,2)*fk2 + beta(4,3)*fk3 + beta(4,4)*fk4);
end

% TO ADD: startup check for ABn 