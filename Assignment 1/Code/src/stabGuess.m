function []=stabGuess(F,degStart,hVec)
%
%
%
%
%
%     

radStart = deg2rad(degStart);
A = [0, 1; -1, 2*cos(radStart)];
S = @(h) max(abs(eig(F(h,A)))) - 1;

fplot(S,hVec, Color = 'k');      % usare per initial guess
hold on
yline(0, LineStyle = '--')
ylim([-1 5])
grid on;    axis padded;    box on;
xlabel('$h$');      ylabel('$S(h)$')      % where S is the stability Slem function
tstring = sprintf('Stability function in alpha = %.0f',degStart);
title (tstring);
drawnow

end