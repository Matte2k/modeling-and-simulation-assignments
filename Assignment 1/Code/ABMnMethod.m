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

%%% ABM3
%   implement generic ABM method in order to have the possiblity to use 2
%   different startup

alpha =  [ 1   2   12  24];

betaP =  [ 1   0   0   0; ... %to check
           3  -1   0   0; ... %to check
           23 -16  5   0; ...
           55 -59  37 -9];    %to check

betaC =  [ 1   0   0   0; ... %to check
           1   1   0   0; ... %to check
           5   8   -1  0; ...
           9   19  -5  1];    %to check


%ABM1
for i = 1 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    xp = x(:,i) + h/alpha(1) * (betaP(1,1)*fk1);
    
    tp = t(i) + h;
    x(:,i+1) = x(:,i) + h/alpha(1) * (betaC(1,1)*f(xp,tp));   % k+1
end


%ABM2
for i = 2 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    xp = x(:,i) + h/alpha(2) * (betaP(2,1)*fk1 + betaP(2,2)*fk2);

    tp = t(i) + h;
    x(:,i+1) = x(:,i) + h/alpha(2) * (betaC(2,1)*f(xp,tp) + betaC(2,2)*fk1);   % k+1
end


%ABM3
for i = 3 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    xp = x(:,i) + h/alpha(3) * (betaP(3,1)*fk1 + betaP(3,2)*fk2 + betaP(3,3)*fk3);
    
    tp = t(i) + h;
    x(:,i+1) = x(:,i) + h/alpha(3) * (betaC(3,1)*f(xp,tp) + betaC(3,2)*fk1 + betaC(3,3)*fk2);   % k+1
end


%ABM4
for i = 4 : (length(t)-1)         % main loop of the method
    fk1 = f(x(:,i)  ,t(i));
    fk2 = f(x(:,i-1),t(i-1));
    fk3 = f(x(:,i-2),t(i-2));
    fk4 = f(x(:,i-3),t(i-3));
    xp = x(:,i) + h/alpha(4) * (betaP(4,1)*fk1 + betaP(4,2)*fk2 + betaP(4,3)*fk3 + betaP(4,4)*fk4);

    tp = t(i) + h;
    x(:,i+1) = x(:,i) + h/alpha(4) * (betaC(4,1)*f(xp,tp) + betaC(4,2)*fk1 + betaC(4,3)*fk2 + betaC(4,4)*fk3);   % k+1
end

% TO ADD: startup check for ABMn 