function [] = stabGuess(F,degStart,hVec)
%STAB GUESS - Plot the Stability problem function in the starting alpha
%
%   Syntax:
%       [] = stabGuess(F,degStart,hVec)
%
%   Input:
%       F,   function(h,A):  corresponding method operator
%       degStart,   double:  starting point of the stability analysis
%       hVec,  double[1,2]:  starting and ending h used in fplot
%
%     

radStart = deg2rad(degStart);
A = [0, 1; -1, 2*cos(radStart)];
S = @(h) max(abs(eig(F(h,A)))) - 1;     % Stability problem definition

fplot(S,hVec, Color = 'k');             % usare per initial guess
hold on;    grid on;    axis padded;    box on; 
ylim(hVec);       yline(0, LineStyle = '--')
xlabel('$h$');    ylabel('$S(h)$')      % where S is the stability function
tstring = sprintf('Stability function in alpha = %.0f',degStart);
title (tstring);
drawnow

end