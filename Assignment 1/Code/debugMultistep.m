clearvars; close all;

f = @(x,t) [-5/2.*(1+8*sin(t)).*x(1); (1-x(1)).*x(2) + x(1)];
x0 = [1 1]';
tmax = 3;
h = 0.01;

figure
hold on

% abOptions2 = abSettings(3,[],[],[]);
% abOptions2.startup = 'RK4';
% [x2,t2,info2] = adamsBashforth(f,x0,tmax,h,abOptions2);

abOptions1 = abSettings(3,[],[],[]);
abOptions1.startup = 'IEX4';
[x1,t1,info1] = adamsBashforth(f,x0,tmax,h,abOptions1);

% abOptions3 = abSettings(3,[],[],[]);
% abOptions3.startup = 'AB';
% [x3,t3,info3] = adamsBashforth(f,x0,tmax,h,abOptions3);

