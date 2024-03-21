clearvars; close all; clc

f = @(x,t) [-x(1); -100*(x(2)-sin(t))+cos(t)];    % function handle of RHS of the ODE
tmax = 1;           % start and ending time
x0 = [1 2]';             % initial guess
theta = 0.5;        % theta param. for the method
h = 0.01;            % time step

options = optimoptions ( 'fsolve', 'Display', 'off' );

m = length ( x0 );              % sys dimension
t = 0:h:tmax;
x = zeros ( m , length(t) );    % solution vector
x(:,1) = x0(:);                 % initial value definition

for i = 1 : (length(t)-1)         % main loop of the method
    tk = t(i);
    xk = x(:,i);

    tinter = tk + theta * h;                    % t_k+1/2
    xinter = xk + theta * h * (f(xk,tk));       % y_k+1/2
    
    tn = tk + h;
    thetaFunction = @(xn) xn - xk - h * theta * f(xk,tk) - (1-theta) * h * f(xn,tn) ;
    
    % DEBUG function plot
    %fplot(theta_residual,[0 2])
    
    % works in both scalar and vectorial cases
    xn = fsolve(thetaFunction, xinter, options);
    
    % works only in scalar case
    %xn = fzero(thetaFunction, xinter);

    % works only in scalar case for the moment
    %xn = newton(thetaFunction, xinter, 'f', 1e-12, 2e3);   
 
    x(:,i+1) = xn;      % solution value in solution vector
end

plot(x(1,:),x(2,:),'bo');

% analytical solution
hold on
plot(exp(-t),2*exp(-100*t)+sin(t),'k');

