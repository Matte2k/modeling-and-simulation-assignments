function[]=stabGuess(F,degStart,hMax)    
    Afun = @(a) [0, 1; -1, 2*cos(a)];
    A = Afun(degStart);
    prob = @(h) max(abs(eig(F(h,A)))) - 1;
    plot([-1 5],[0 0],'k--')
    hold on
    fplot(prob,[0 hMax]);      % usare per initial guess
    ylim([-1 5])
    drawnow
end