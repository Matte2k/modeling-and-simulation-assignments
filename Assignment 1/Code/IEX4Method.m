clearvars; close all; clc

addpath(genpath('src\'));
%%%%                        --- TO FIX ---

% input system 
b11 = -180.5;     b12 = 219.5;
b21 = 179.5;      b22 = -220.5;
B = [b11 b12; b21 b22];
I = eye(2);

alpha = [-1/6, 4, -27/2, 32/3];

f = @(x,t) [b11.*x(1) + b12.*x(2); b21.*x(1) + b22.*x(2)];    % function handle of RHS of the ODE
tmax = 1;           % start and ending time
x0 = [1 1]';        % initial guess
h = 0.1;            % time step

options = optimoptions ( 'fsolve', 'Display', 'final' );

m = length ( x0 );              % sys dimension
t = 0:h:tmax;
x = zeros ( m , length(t) );    % solution vector
x(:,1) = x0(:);                 % initial value definition

matInput = false;

for i = 1 : (length(t)-1)         % main loop of the method
    
    if matInput == false
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


    elseif matInput == true
        xp1 = (I + h*B) * x(:,i);
        xp2 = (I + h*B).^2 * x(:,i);
        xp3 = (I + h*B).^3 * x(:,i);
        xp4 = (I + h*B).^4 * x(:,i);
  
    end
    
    x(:,i+1) = alpha(1) * xp1 + alpha(2) * xp2 + alpha(3) * xp3 + alpha(4) * xp4;
end

figure
plot(x(1,:),x(2,:),'bo');
hold on

% analytical solution
%figure('Name','anal')
xAnal1 = @(t) exp(B(1,:).*t) * [1;1];
xAnal2 = @(t) exp(B(2,:).*t) * [1;1];

xAn1Point(1) = 1;
xAn2Point(1) = 1;
for i = 2 : (length(t))
    tNow = t(i);
    xAn1Point(i) = xAnal1(tNow);
    xAn2Point(i) = xAnal2(tNow);
end

plot(xAn1Point,xAn2Point,'xr')

%plot(xAnal1,xAnal2,'k');
% xlim([-0.5 0.5])
% ylim([-1 100])

%figure('Name','rk4')
[x_rk,t_rk,info] = rk4(f,x0 ,tmax ,0.1);
plot(x_rk(1,:),x_rk(2,:))
