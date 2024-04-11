function [xF,yF,hvec,stabEdge]=stabRegion(F,stepNum,degVec,guess,visualConfig)
%
%
%
%
%
%
%
%
%

if nargin < 5
    visualConfig = true;
end

if length(guess) > 2
    error('guess dimension incorrect. Insert a valid guess')
end


%%% Variables initialization
Afun = @(a) [0, 1; -1, 2*cos(a)];
alphaLim = deg2rad(degVec);
alphaVec = linspace(alphaLim(1),alphaLim(2),stepNum);
hvec = zeros(1,length(alphaVec));
xF   = zeros(1,length(alphaVec));
yF   = zeros(1,length(alphaVec));



for i=1:length(alphaVec)
    
    alpha = alphaVec(i);
    A = Afun(alpha);
    S = @(h) max(abs(eig(F(h,A)))) - 1;    
    [hvec(i),~,conv]=fzero(S, guess);
    
    if length(guess) == 1
        guess = hvec(i);        
    elseif length(guess) == 2
        guess(2) = hvec(i);       
    end

    % Eigenvalue of discrete Slem
    eig_iterVec = eig(A);
    xF(:,i) = hvec(i)*real(eig_iterVec(1));
    yF(:,i) = hvec(i)*imag(eig_iterVec(1));
    
    if conv < 1
        warning('Implicit equation at the step %d has not been solved correctly',i)
    end
end

% stability region coordinate vectors
stabEdge = [xF,  flip(xF,2);yF, -flip(yF,2)];


if visualConfig == true
    axis equal;     grid on;    box on;     hold on
    ax = gca;
        ax.Layer = 'top';
    % final stability region computed plot
    plot(stabEdge(1,:), stabEdge(2,:),'k')
    if degVec(1) ~= 0
        fill(stabEdge(1,:), stabEdge(2,:),[0.80 0.80 0.80])
    elseif degVec(1) == 0
        ax.Color = [0.80 0.80 0.80];
        fill([-inf -inf inf inf],[-inf inf inf -inf],[0.80 0.80 0.80])
        fill(stabEdge(1,:), stabEdge(2,:),'w')
        grid on;
    end
    xline(0,'--');                  yline(0,'--')
    xlabel('$Re\{h\lambda\}$');     ylabel('$Im\{h\lambda\}$');
    title ('Stability region')

end


end
