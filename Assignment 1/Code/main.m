% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assugnment # 1
% Author: Matteo Baio 10667431

addpath(genpath('src\'));

%% Ex 1 - CELL
clearvars; close all; clc;
%clearvars; close all; 

% nota il sistema dev'essere quadrato (n variabili ed n equazioni)

%%% FUNCTION INPUT
f = @(x1,x2) [x2.^2 + x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
xGuess = [-1 1]';
toll = 1e-10;
nmax = 200;
feval = 0;
dimSys = 2;  % DA CAPIRE SE POSSO EVITARE DI METTERLO COME INPUT

% pick up a method:
method = 's';

% set flag
printFlag = true;
plotFlag = true;


%%% FUNCTION BODY  - to add standard values for nmax, toll and initial guess if not specified

%%% Initialization
timerStart = tic;   % timer start
err = toll + 1;     % error
k = 0;              % iterator
errVect = [];       % error vector
xkVect = [];        % solutions vector
xGuessCell = num2cell(xGuess);  % inital guess cell array


%%% Newton method body
% Check if initial guess is already a solution
if f(xGuessCell{:,1}) < toll*ones(length(xGuess),1)
    solution = xGuess;
    converge = true;
    elapsedTime = toc(timerStart);

    info.iteration = 0;
    info.solutionVector = solution;
    info.errorVector = f(xGuessCell{:,1});
    info.timeCost = elapsedTime;
    info.fevalCost = dimSys;
    info.methodUsed = 'none';

    if printFlag == true
        fprintf('Function zeros are equal to the initial guess\n');
        fprintf('The error %.3d is less than the tollerance %.3d \n', ...
                   errVect(end), toll);

        solPrint = sprintf(' %.4f\n',solution);
        fprintf('\n Zeros computed:\n%s \n',solPrint)
    end
    
    % function in this "lucky" case ends here
    return

else
    
    feval = feval + dimSys;     % function evaluation used to check if the guess is already a solution
    xk = xGuess;                % start point of the iteration process
   
    % symbolic manipulation to compute jacobian matrix
    if method == 's'
        fSym = sym(f);
        JSym = jacobian(fSym);
        J = matlabFunction(JSym);
    end

    % Iterative loop to find the solution starts
    while k < nmax && err > toll

        if method == 's'            %  - TO SUBSTITUTE IT WITH A SWITCH CASE - 
            % function eval
            xkCell = num2cell(xk);
            fk = f(xkCell{:,1});          % feval + z
            feval = feval + dimSys;                         % TO CONTROL

            Jk = J(xkCell{:,1});          % feval + z^z
            feval = feval + dimSys^(dimSys);                % TO CONTROL

            methodString = 'SYMBOLIC MATH';  % method used saved for the output

        elseif method == 'f'
            epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));     % dimSys,1 vector
            Jk = zeros(dimSys);

            xkCell = num2cell(xk);
            fk = f(xkCell{:,1});
            feval = feval + dimSys;                         % TO CONTROL

            for i=1:dimSys
                epsilonVect = (zeros(dimSys,1));
                epsilonVect(i,1) = epsilon(i,1);
                
                xkForward =  xk + epsilonVect;
                xkForwardCell = num2cell(xkForward);
                fkForward = f(xkForwardCell{:,1});
                feval = feval + dimSys;                     % TO CONTROL

                Jk(:,i) = (fkForward - fk) / epsilon(i,1); 
            end

            methodString = 'FORWARD DIFF';  % method used saved for the output

        elseif method == 'c'
            epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));
            Jk = zeros(dimSys);

            xkCell = num2cell(xk);
            fk = f(xkCell{:,1});
            feval = feval + dimSys;                         % TO CONTROL

            for i=1:dimSys
                epsilonVect = (zeros(dimSys,1));
                epsilonVect(i,1) = epsilon(i,1);
                
                xkForward =  xk + epsilonVect;
                xkForwardCell = num2cell(xkForward);
                fkForward = f(xkForwardCell{:,1});
                feval = feval + dimSys;                     % TO CONTROL

                xkBackward =  xk - epsilonVect;
                xkBackwardCell = num2cell(xkBackward);
                fkBackward = f(xkBackwardCell{:,1});
                feval = feval + dimSys;                     % TO CONTROL
                
                Jk(:,i) = (fkForward - fkBackward) / (2*epsilon(i,1)); 
            end        
 
            methodString = 'CENTERED DIFF';  % method used saved for the output

        end

        % check if f is sufficiently well-behaved
        if rank(Jk) == zeros(length(xGuess))
            error('Jk has become zero, f not sufficiently well-behaved near solution\n');

        else
            y = Jk \ (-fk);
            xNew = xk + y;              % solutions for the k iteration
            err = norm(y);              % errore commesso differenza iterate

            errVect = [errVect, err];   % error vector update
            xkVect = [xkVect, xNew];    % solution vector update
            k = k + 1;                  % iteration index update
            xk = xNew;                  % starting point update
        end
    end
end

elapsedTime = toc(timerStart);   % timer stop

%%% Output organization
solution = xkVect(:,end);   % final result
if err > toll           % convergence flag 
    converge = false;
else
    converge = true;
end

info = struct;      % info struct for detailed output
    info.iteration = k;
    info.solutionVector = xkVect;
    info.errorVector = errVect;
    info.timeCost = elapsedTime;
    info.fevalCost = feval;
    info.methodUsed = methodString;

% main output print on command window
if printFlag == true
    if converge == true
        fprintf('Function zeros has been successfully computed using %s method \n',methodString);
        fprintf('After %d interations the error %.3d is less than the tollerance %.3d \n', ...
                    k, errVect(end), toll);
    else
        fprintf('Function zeros has NOT been computed using %s method \n',methodString);
        fprintf('After %d interations the error %.3d is still grater than the tollerance %.3d \n', ...
                    k, errVect(end), toll);
    end
    solPrint = sprintf(' %.4f\n',solution);
    fprintf('\n Zeros computed:\n%s \n',solPrint)
end

% convergence plot 
if plotFlag == true
    %to write plot
end


%%% DEBUG - compare with built in function
% disp('----------------------------------------')
% fMatlab = @(x) [x(2).^2 + x(1) - 2; -x(1).^2 + x(2) + 10];
% fsolve(fMatlab, xGuess)

%% prova

f = @(x1,x2) [x2.^2 + x1 - 2; -x1.^2 + x2 + 10];    % # di funzioni = z
xGuess = [-1 1]';
toll = 1e-10;

% figure()
% grid on
% hold on

[solution,converge,info] = newton(f, xGuess, 'f', toll, 200, false, true);
% [solution,converge,info] = newton(f, xGuess, toll, 'f', 200)
% [solution,converge,info] = newton(f, xGuess, toll, 'c', 200)



%% Ex 1 - MODEL WITHOUT CELL ARRAY (incomplete)
% clearvars; close all; clc;
% 
% % data input
% f = @(x) [x(2).^2 + x(1) - 2; -x(1).^2 + x(2) + 10];
% xGuess = [-1 1]';
% toll = 1e-10;
% 
% % pick up a method:
% method = 's';
% 
% %%% FUNCTION BODY
% err = toll + 1;     % inizializzazione errore
% k = 0;              % inizializzazione iteratore
% xkVect = [];         % inizializzazione vettore soluzioni
% xk = xGuess;
% nmax = 200;
% 
% if method == 's'
%     % symbolic manipulation
%     syms x1 x2
%     x = [x1;x2];
%     fSym = f(x); 
%     JSym = jacobian(fSym);
%     J = matlabFunction(JSym,'Vars',{[x1; x2]});
% end
% 
% while k < nmax && err > toll
% 
%     if method == 's'      
%         % function eval
%         fk = f(xk);
%         Jk = J(xk);
% 
%     elseif method == 'f'
%         epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));
% 
%     elseif method == 'c'
%         epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));
% 
%     end
% 
% 
% 
%     if fk == zeros(length(xGuess),1)
%         disp('capiamo...')  
%     elseif Jk == 0
%         error('Arresto per azzeramento di dfun');
%     else
%         y = Jk \ (-fk);
%         xNew = xk + y;             % calcolo soluzione k-esima iterata
%         err = norm(y);           % errore commesso differenza iterate
%         xkVect = [xkVect, xNew];     % aggiornamento soluzione
%         k = k + 1;                          
%         xk = xNew;                % aggiornamento k-esima iterata
%     end
% 
%    % AGGIUNGI errore QUANDO NON CONVERGE
% 
% end
% 
% %%% FUNCTION OUTPUT
% zero = xkVect(:,end)
% xkVect;
% iter=k;
% error = 0; % to be defined error vector
% 
% %verifica
% fMatlab = @(x) [x(2).^2 + x(1) - 2; -x(1).^2 + x(2) + 10];
% fsolve(fMatlab, xGuess)
