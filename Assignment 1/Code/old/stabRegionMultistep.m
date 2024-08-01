clearvars; close all; clc;
graphicSettings;

dimSys = 2;
stepNum = 100;

%%% RK2 STABILITY PROBLEM
figure('Name','AB2');
    t1 = tiledlayout(1,2);
    tileTitle = title(t1,'AB2 stability analysis');
    tileTitle.Interpreter = 'latex';

nexttile
    [AB2.F,AB2.guess,AB2.degVec,AB2.hVec]=selectOp('AB2',dimSys);
    stabGuess(AB2.F,AB2.degVec(1),AB2.hVec);

nexttile
    [AB2.xF,AB2.yF,AB2.hvec,AB2.stabEdge]=stabRegion(AB2.F,stepNum,AB2.degVec,AB2.guess);
    xlim([-4 1])
    
leg = legend('$AB2$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
% -add graphic export here-

fprintf('Solution of the problem in alpha = pi using RK2 is:\n h = %.4f \n\n',AB2.hvec(1));

%%


%   F_m1 = (I-h*A + ((1-params)*h*A)^2/2) \ (I);


%%
clearvars; close all; clc;
graphicSettings;

dimSys = 2;
stepNum = 100;

%%% RK2 STABILITY PROBLEM
figure('Name','AB2');
    t1 = tiledlayout(1,2);
    tileTitle = title(t1,'AB2 stability analysis');
    tileTitle.Interpreter = 'latex';

nexttile
    [AB2.F,AB2.guess,AB2.degVec,AB2.hVec]=selectOp('BI2',dimSys,0);
    stabGuess(AB2.F,AB2.degVec(1),AB2.hVec);

nexttile
    [AB2.xF,AB2.yF,AB2.hvec,AB2.stabEdge]=stabRegion(AB2.F,stepNum,AB2.degVec,AB2.guess);
    xlim([-4 1])
    
leg = legend('$AB2$','$Stable$','Orientation', 'Horizontal');
    leg.Layout.Tile = 'south';
% -add graphic export here-

fprintf('Solution of the problem in alpha = pi using RK2 is:\n h = %.4f \n\n',AB2.hvec(1));