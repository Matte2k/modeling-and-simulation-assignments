function [Z1mesh,Z2mesh,t] = zerosGuess(f1,f2,meshMin,meshMax,meshStep,visualConfig)
%ZEROS GUESS - plot the two functions on a defined mesh to find an educated
% guess for 2x2 Newton method
%
%   Syntax:
%       [Z1mesh,Z2mesh,t] = zerosGuess(f1,f2,meshMin,meshMax,visualConfig)
%
%   Input:
%       f1,  function(x1,...,xn):  function 1 for which to compute zeros
%       f2,  function(x1,...,xn):  function 2 for which to compute zeros
%       meshMin(*),       double:  mesh lower limit
%       meshMax(*),       double:  mesh upper limit
%       meshStep(*),      double:  mesh step between to values
%       visualConfig(*),    bool:  set visual output of the function
%
%   Output:
%       Z1mesh,   double[n,m]:  mesh evaluated on f1
%       Z2mesh,   double[n,m]:  mesh evaluated on f2
%       t,   TiledChartLayout:  visual output object
%
%   Default settings for optional input (*):
%       meshMin:       set equal to -5 by default
%       meshMax:       set equal to -5 by default
%       meshStep:      set equal to 0.1 by defaylt
%       visualConfig:  set as true by default
%


    %%% Default value for optional input
    if nargin < 6
        visualConfig = true;
        if nargin < 5
            meshStep = 0.1;
            if nargin < 4
                meshMax = 5;
                if nargin < 3
                    meshMin = -5;
                end
            end
        end
    end
    
    
    %%% mesh definition and evaluation
    [X1mesh,X2mesh] = meshgrid(meshMin:meshStep:meshMax);
    Z1mesh = f1(X1mesh,X2mesh);
    Z2mesh = f2(X1mesh,X2mesh);
    
    
    if visualConfig == true
        fig = figure();
            fig.Name = 'Function plot';
            %fig.Position = [1 1 4000 2000];    % to find values
        t = tiledlayout(1,2);
        tString = title(t,'Initial guess graphic search');
        tString.Interpreter = 'latex';
        
        %%% plot 3d of the surface
        nexttile
        hold on;    box  on;    grid on;    axis padded
        view(3)
        
        % surface plot
        s0 = surf(X1mesh,X2mesh,zeros(length(Z1mesh)));
            s0.LineStyle = "none";
            s0.FaceColor = 'k';
            s0.FaceAlpha = 0.2;
            
        s1 = surf(X1mesh,X2mesh,Z1mesh);
            s1.LineStyle = "none";
            s1.FaceColor = 'r';
            s1.FaceAlpha = 0.4;
        
        s2 = surf(X1mesh,X2mesh,Z2mesh);
            s2.LineStyle = "none";    
            s2.FaceColor = 'b';
            s2.FaceAlpha = 0.4;
        
        % % versors
        % x1ax = quiver3(0,0,0,1,0,0);
        %     x1ax.LineWidth = 2;
        %     x1ax.Color = 'k';
        % x2ax = quiver3(0,0,0,0,1,0);
        %     x2ax.LineWidth = 2;
        %     x2ax.Color = 'k';
        % zax = quiver3(0,0,0,0,0,1);
        %     zax.LineWidth = 2;
        %     zax.Color = 'k';
        
        % labels
        xlabel('$x_1$');   ylabel('$x_2$');
        title('Surface plot')
        legend('$zero$','$f_1$','$f_2$')
        
        
        %%% contour plot of the zero level
        nexttile
        hold on;    box  on;    grid on;    axis padded
        
        % contour plot level 0
        [~,c1] = contour(X1mesh,X2mesh,Z1mesh,[0 0]);
            c1.EdgeColor = 'r';
        
        [~,c2] = contour(X1mesh,X2mesh,Z2mesh,[0 0]);
            c2.EdgeColor = 'b';
        
        % labels
        xlabel('$x_1$');   ylabel('$x_2$');  zlabel('$f(x_1,x_2)$')
        title('Contour plot')
        legend('$f_1=0$','$f_2=0$')
    
    end

end