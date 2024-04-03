clearvars; close all;
addpath(genpath('src\'));

% Data
f = @(x,t) [-5/2.*(1+8*sin(t)).*x(1); (1-x(1)).*x(2) + x(1)];
x0 = [1 1]';
tmax = 1;
h = 0.1;


%% AB
figure("Name",'Adams Bashforth')
hold on
grid on

abOptionsSS = abSettings(3,[],[],[]);
abOptionsSS.startup = 'RK4';
[ab.ss.x,ab.ss.t,ab.ss.info] = adamsBashforth(f,x0,tmax,h,abOptionsSS);

abOptionsMS = abSettings(3,[],[],[]);
abOptionsMS.startup = 'AB';
[ab.ms.x,ab.ms.t,ab.ms.info] = adamsBashforth(f,x0,tmax,h,abOptionsMS);


%% AM
figure("Name",'Adams Moulton')
hold on
grid on

amOptionsSS = amSettings(3,[],[],[],[]);
amOptionsSS.startup = 'RK4';
[am.ss.x,am.ss.t,am.ss.info] = adamsMoulton(f,x0,tmax,h,amOptionsSS);

amOptionsMS = amSettings(3,[],[],[],[]);
amOptionsMS.startup = 'AM';
[am.ms.x,am.ms.t,am.ms.info] = adamsMoulton(f,x0,tmax,h,amOptionsMS);


%% ABM
figure("Name",'Adams Bashforth Moulton')
hold on
grid on

abmOptionsSS = abmSettings(3,[],[],[],[]);
abmOptionsSS.startup = 'RK4';
[abm.ss.x,abm.ss.t,abm.ss.info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptionsSS);

abmOptionsMS = abmSettings(3,[],[],[],[]);
abmOptionsMS.startup = 'ABM';
[abm.ms.x,abm.ms.t,abm.ms.info] = adamsBashforthMoulton(f,x0,tmax,h,abmOptionsMS);


%% BDF
figure("Name",'Backward Difference Formula')
hold on
grid on

bdfOptionsSS = bdfSettings(3,[],[],[],[],[]);
bdfOptionsSS.startup = 'RK4';
[bdf.ss.x,bdf.ss.t,bdf.ss.info] = backwardDifferenceFormula(f,x0,tmax,h,bdfOptionsSS);

bdfOptionsMS = bdfSettings(3,[],[],[],[],[]);
bdfOptionsMS.startup = 'BDF';
[bdf.ms.x,bdf.ms.t,bdf.ms.info] = backwardDifferenceFormula(f,x0,tmax,h,bdfOptionsMS);


%% COMPARE
figure("Name",'Backward Difference Formula')
hold on
grid on
plot(am.ss.t,am.ss.x,'k-')
plot(ab.ss.t,ab.ss.x,'r-')
plot(abm.ss.t,abm.ss.x,'g-')
plot(bdf.ss.t,bdf.ss.x,'b-')
legend('AM','AM','AB','AB','ABM','ABM','BDF','BDF','location','best')