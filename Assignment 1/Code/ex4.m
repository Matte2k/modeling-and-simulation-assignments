clearvars;  close all;   clc

Afun = @(a) [0, 1; -1, 2*cos(a)];
x0 = [1 1]';
tmax = 1;
tollVec = [1e-3 1e-4 1e-5 1e-6];         % output di altra function
color = 'rgbm'; 

%%% INPUT
degStart = 180;
degEnd = 0;
stepNum = 100;
method = 'fzero';   % default
dimSys = 2;
[RKcell{1}]=selectOp('RK1',dimSys);
[RKcell{2}]=selectOp('RK2',dimSys);
[RKcell{3}]=selectOp('RK4',dimSys);


%%% Variables initialization
alphalim = deg2rad([degStart degEnd]);
alphaVec = linspace(alphalim(1),alphalim(2),stepNum);

%fsolve options settings
hMax = 1;

% Guess finding
for ord = 1:length(RKcell)
    figure("Name",'Initial Guess')
    axis equal;     grid on;    hold on;
    RK = RKcell{ord};
    for tols = 1:length(tollVec)

        alphaStart = alphaVec(1);
        A = Afun(alphaStart);
        analSol = expm(A*tmax)*x0;

        nstep = @(h) (tmax-0)/h;
        rkSol = @(h) (RK(h,A))^nstep(h) * x0;
        prob = @(h) norm((analSol-rkSol(h)),inf)-tollVec(tols);

        % to check mark on graphic prob(toll) and prob(1) 

        plot([-1 hMax],[0 0],'k--')
        hold on
        fplot(prob,[0 hMax]);      % MODO PIU' INTELLIGENTE??
        ylim([-1 hMax])
        drawnow

    end
end


for ord = 1:length(RKcell)
    figure
    axis equal;     grid on;    hold on
    RK = RKcell{ord};
    
    for tols=1:length(tollVec)
        
        toll = tollVec(tols);
        guess = [toll 1];     % INPUT
    
        for i=1:length(alphaVec)
            
            alpha = alphaVec(i);
            A = Afun(alpha);
    
            analSol = expm(A*tmax)*x0;
            
            nstep = @(h) (tmax-0)/h;
            rkSol = @(h) (RK(h,A))^nstep(h) * x0;
    
            prob = @(h) norm((analSol-rkSol(h)),inf)-toll;    
            
            %[hvec(i),fval(i),conv(i)]=fsolve(prob,guess);
            [hvec(i),fval(i),conv(i)]=fzero(prob,guess);
            guess = hvec(i);    % initial guess update
                   
            % Eigenvalue of discrete problem
            eig_iterVec = eig(A);
            xF(i) = hvec(i)*real(eig_iterVec(1));
            yF(i) = hvec(i)*imag(eig_iterVec(1));        
        end
    
        % final stability region computed plot
        plot(xF,yF,'color',color(tols))
        fill(xF,yF,color(tols),'FaceAlpha', 0.3)
        plot(xF,-yF,'color',color(tols))
        fill(xF,-yF,color(tols),'FaceAlpha', 0.3)
    end
end
