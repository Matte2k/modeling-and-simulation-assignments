clearvars; close all; clc

addpath(genpath('src\'));

% input system 
b11 = -180.5;     b12 = 219.5;
b21 = 179.5;      b22 = -220.5;
B = [b11 b12; b21 b22];
I = eye(2);

alpha = [-1/6, 4, -27/2, 32/3];

f = @(x,t) [b11.*x(1) + b12.*x(2); b21.*x(1) + b22.*x(2)];    % function handle of RHS of the ODE
tmax = 1;           % start and ending time
x0 = [1 1]';        % initial guess
h = 0.01;            % time step

options = optimoptions ( 'fsolve', 'Display', 'final' );

m = length ( x0 );              % sys dimension
t = 0:h:tmax;
x = zeros ( m , length(t) );    % solution vector
x(:,1) = x0(:);                 % initial value definition

figure(1)
grid on
xlabel('x1')
ylabel('x2')
hold on

for i = 1 : (length(t)-1)         % main loop of the method

    %1st predictor
    tk1 = t(i) + h;
    fk1 = @(k1) x(:,i) + h * f(k1,tk1) - k1;
    xp1 = fsolve(fk1, x(:,i), options);         % initial guess xk

    %2nd predictor
    tk2a = t(i) + h/2;
    fk2a = @(k2a) x(:,i) + h/2 * f(k2a,tk2a) - k2a;
    k2a = fsolve(fk2a, x(:,i), options);        % initial guess xk

    tk2 = tk1;
    fk2 = @(k2) k2a + h/2 * f(k2,tk2) - k2;
    xp2 = fsolve(fk2, k2a, options);            % initial guess k2a

    %3rd predictor
    tk3a = t(i) + h/3;
    fk3a = @(k3a) x(:,i) + h/3 * f(k3a,tk3a) - k3a;
    k3a = fsolve(fk3a, x(:,i), options);        % initial guess xk

    tk3b = t(i) + h * 2/3;
    fk3b = @(k3b) k3a + h/3 * f(k3b,tk3b) - k3b;
    k3b = fsolve(fk3b, k3a, options);           % initial guess k3a

    tk3 = tk1;
    fk3 = @(k3) k3b + h/3 * f(k3,tk3) - k3;
    xp3 = fsolve(fk3, k3b, options);            % initial guess k3b

    %4th predictor
    tk4a = t(i) + h/4;
    fk4a = @(k4a) x(:,i) + h/4 * f(k4a,tk4a) - k4a;
    k4a = fsolve(fk4a, x(:,i), options);        % initial guess xk

    tk4b = t(i) + h * 2/4;
    fk4b = @(k4b) k4a + h/4 * f(k4b,tk4b) - k4b;
    k4b = fsolve(fk4b, k4a, options);           % initial guess k4a

    tk4c = t(i) + h * 3/4;
    fk4c = @(k4c) k4b + h/4 * f(k4c,tk4c) - k4c;
    k4c = fsolve(fk4c, k4b, options);           % initial guess k4b

    tk4 = tk1;
    fk4 = @(k4) k4c + h/4 * f(k4,tk4) - k4;
    xp4 = fsolve(fk4, k4c, options);            % initial guess k4c

    %corrector
    x(:,i+1) = alpha(1) * xp1 + alpha(2) * xp2 + alpha(3) * xp3 + alpha(4) * xp4;

end

plot(x(1,:),x(2,:),'kx');

% analytical solution
xAnal = @(t) expm(B.*t) * [1;1];    % 'expm' since is a matrix exponential
xAnalPoint(:,1) = [1;1];
for i = 2 : (length(t))
    tNow = t(i);
    xAnalPoint(:,i) = xAnal(tNow);
end
plot(xAnalPoint(1,:),xAnalPoint(2,:),'r')


% ODE verification
fODE = @(t,x) [b11.*x(1) + b12.*x(2); b21.*x(1) + b22.*x(2)];    % function handle of RHS of the ODE
[t,y] = ode15s(fODE,[0 1],x0);
plot(y(:,1),y(:,2),'bo')


% richardson and RK matrix method 
x_mat(:,1) = x0(:);
x_matRK(:,1) = x0(:); 
for i = 1 : (length(t)-1)         % main loop of the method
    xp1_mat = (I + h*B) * x(:,i);
    xp2_mat = (I + h*B).^2 * x(:,i);
    xp3_mat = (I + h*B).^3 * x(:,i);
    xp4_mat = (I + h*B).^4 * x(:,i);

    x_mat(:,i+1) = alpha(1) * xp1_mat + alpha(2) * xp2_mat + alpha(3) * xp3_mat + alpha(4) * xp4_mat;
    x_matRK(:,i+1) = (I + (h*B) + 1/2*(h*B).^2 + 1/6*(h*B).^3 + 1/24*(h*B).^4) * x(:,i);
end
% plot(x_mat(1,:),x_mat(2,:),'go');
% plot(x_matRK(1,:),x_matRK(2,:),'mo');

legend('IEX4','anal','ode','richard','rk')

