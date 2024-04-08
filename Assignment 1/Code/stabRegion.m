% Stability region computation example
clearvars
close all
clc

addpath(genpath('src\'));

% % RK4 and RK2, FUNCTION CHE RACCOGLIE OPERATORI E DA CUI SELEZIONI
% I = eye(2);
Afun = @(a) [0, 1; -1, 2*cos(a)];
% F4 = @(h,A) I + (h*A) + 0.5*(h*A)^2 + 1/6*(h*A)^3 + 1/24*(h*A)^4;
% F3 = @(h,A) I + (h*A) + 1/2*(h*A)^2 + 1/6*(h*A)^3;    
% F2 = @(h,A) I + (h*A) + 0.5*(h*A)^2;
% F1 = @(h,A) I + (h*A);
% theta = 0.4;
% BI2 = @(h,A) (I-(1-theta)*h*A+((1-theta)*h*A)^2/2) \ (I+theta*h*A+(theta*h*A)^2/2);
% IEX4
% esercizio 4

%%% F definition
F = cell(1,4);
dimSys = 2;
[F{1}]=selectOp('BI2',dimSys,0.7);
[F{2}]=selectOp('BI2',dimSys,0.9);
[F{3}]=selectOp('RK3',dimSys,0.3);
[F{4}]=selectOp('RK4',dimSys,0.2);

%%% INPUT
degStart = 180;
degEnd = 0;
stepNum = 50;
method = 'fzero';   % default

%%% Variables initialization
color = ['b' 'r' 'g' 'm'];
alphalim = deg2rad([degStart degEnd]);
alphaVec = linspace(alphalim(1),alphalim(2),stepNum);


%fsolve options settings
options = optimset('Display','iter');
hMax = 15;
figure("Name",'Initial Guess')
axis equal;     grid on;    hold on
for ord=1:2
    guessFinder(F{ord},degStart,hMax);
end
guessVec = [10 10 2 2];         % output di altra function


figure
axis equal;     grid on;    hold on

for ord=1:2
    guess = guessVec(ord);

    for i=1:length(alphaVec)
        Fnow = F{ord};

        alpha = alphaVec(i);
        A = Afun(alpha);
        prob = @(h) max(abs(eig(Fnow(h,A)))) - 1;    
        
        %%% H COMPUTATIONS
        %[hvec(i),conv(i)] = newton(prob, guess , 'f' ,1e-7, 1e3);
        %[hvec(i),fval(i),conv(i)]=fzero(prob,guess,options);
        
        %[hvec(i),fval(i),conv(i)]=fsolve(prob,guess);
        [hvec(i),fval(i),conv(i)]=fzero(prob,guess);
        guess = hvec(i);    % initial guess update

        % % Eigenvalue of continuous problem
        eig_iterVec = eig(A);
        xA(i) = real(eig_iterVec(1));
        yA(i) = imag(eig_iterVec(1));
        % plot(xA(i),yA(i),Marker="o",MarkerEdgeColor='k')
        
        % Eigenvalue of discrete problem
        xF(i) = hvec(i)*real(eig_iterVec(1));
        yF(i) = hvec(i)*imag(eig_iterVec(1));
        %plot(xF(i),yF(i),Marker="+",MarkerEdgeColor='r')    
    end

    % final stability region computed plot
    plot(xF,yF,'color',color(ord))
    fill(xF,yF,color(ord),'FaceAlpha', 0.3)
    plot(xF,-yF,'color',color(ord))
    fill(xF,-yF,color(ord),'FaceAlpha', 0.3)
end


