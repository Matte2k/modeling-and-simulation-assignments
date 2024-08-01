function [xF,yF,hvec,stabEdge] = stabRegion(F,stepNum,degVec,guess,visualConfig)
%STAB REGION - Finds the stability region associated to the operator F(h,A)
%
%    Syntax:
%       [xF,yF,hvec,stabEdge] = stabRegion(F,stepNum,degVec,guess,visualConfig)
%
%   Input:
%       F,       function(h,A):  method operator to compute stability region
%       stepNum,        double:  number of intermidiate step between alpha limit 
%       degVec,    double[1,2]:  starting and ending alpha to compute stability region
%       guess,          double:  eductaed initial guess to find stability region
%       visualConfig(*),  bool:  set as true to plot solutions 
%
%   Output:
%       xF,       double[1, n]:  x coordinates of the h*lambda plane points between degVec values
%       yF,       double[1, n]:  y coordinates of the h*lambda plane points between degVec values
%       hVec,     double[1, n]:  h computed for the stability problem between degVec values
%       stabEdge, double[2,2n]:  stability region contour coordinates in all the h*lambda plane
%
%   Default settings for optional input (*):
%       visualConfig:  set as true by default
%

    
    %%% Optional input definition
    if nargin < 5
        visualConfig = true;
    end
    
    %%% Variables initialization
    Afun = @(a) [0, 1; -1, 2*cos(a)];
    alphaLim = deg2rad(degVec);
    alphaVec = linspace(alphaLim(1),alphaLim(2),stepNum);
    hvec = zeros(1,length(alphaVec));
    xF   = zeros(1,length(alphaVec));
    yF   = zeros(1,length(alphaVec));
    
    %%% Stability analysis loop
    for i=1:length(alphaVec)        
        
        alpha = alphaVec(i);
        A = Afun(alpha);                        % Local A matrix value
        S = @(h) max(abs(eig(F(h,A)))) - 1;     % Stability problem definition
        [hvec(i),~,conv]=fzero(S, guess);       % Stability problem resolution
              
        guess = hvec(i);        % Educated gusess based on the current solution       
        eig_iterVec = eig(A);   % Eigenvalue of continuos system
        xF(:,i) = hvec(i)*real(eig_iterVec(1)); % h*lambda x coordinates
        yF(:,i) = hvec(i)*imag(eig_iterVec(1)); % h*lambda y coordinates
        
        % convergence check
        if conv < 1
            warning('Implicit equation at the step %d has not been solved correctly',i)
        end
    end
    
    % Entire stability region coordinate vectors
    stabEdge = [xF,  flip(xF,2);yF, -flip(yF,2)];
    
    if visualConfig == true
        axis equal;     grid on;    box on;     hold on
        ax = gca;       ax.Layer = 'top';
            
        % final stability region computed plot
        plot(stabEdge(1,:), stabEdge(2,:),'k')
        if degVec(1) ~= 0       % explicit method case
            fill(stabEdge(1,:), stabEdge(2,:),[0.80 0.80 0.80])
        elseif degVec(1) == 0   % implicit method case
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