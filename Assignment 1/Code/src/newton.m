function [solution,converge,info] = newton(f,xGuess,method,toll,nmax,visualConfig)
%NEWTON - Compute zeros of given function using the Newton's method
%
%   Syntax:
%       [solution,converge,info] = newton(f,xGuess,method,toll,nmax,visualConfig)
%
%   Input:
%       f,   function(x1,...,xn):  function for which to compute zeros
%       xGuess,      double[n,1]:  generic initial guess 
%       method(*),          char:  chose between:
%               - 's' compute jacobian using symbolic math
%               - 'f' compute jacobian using forward finite difference
%               - 'c' compute jacobian using centered finite difference
%       toll(*),          double:  tollerance on the error between iterations
%       nmax(*),          double:  maximum number of iterations
%       visualConfig(*),  struct:  set visual output of the function
%               - visualConfig.print: display solutions and some info
%               - visualConfig.plot: plot the convergence plot
%
%   Output:
%       solution,           double:  solution vector
%       convergence,   double[1,m]:  time istant associated to solutions
%       info,               struct:  information on method used:
%           - info.iteration,           double:  number of iteration (m)
%           - info.solutionVector, double[n,m]:  solution vector for each iteration
%           - info.errorVector,    double[1,m]:  error for each iteration
%           - info.timeCost,            double:  time spent
%           - info.fevalCost,           double:  # of function evaluations
%           - info.methodUsed,            char:  string describing the method used
%
%   Default settings for optional input (*):
%       method:        set as 'f' by default
%       toll:          set as 1e-12 by default
%       nmax:          set as 1e3 by default
%       visualConfig:  set as true by default
%


    %%% Prelimirary check
    if nargin < 6                       % default options without visual config flags
        visualConfig.print = true;
        visualConfig.plot = true;
    end

    if nargin < 5 || isempty(nmax)      % default max number of iteration
        nmax = 1e3;
    end

    if nargin < 4 || isempty(toll)      % default tollerance
        toll = 1e-12;
    end

    if method ~= 's' && method ~= 'f' && method ~= 'c'
        error ('Method not valid, please insert a valid method as input')
    end


    %%% Initialization
    timerStart = tic;               % timer start
    feval = 0;                      % function evaluation counter starts
    dimSys = length(xGuess);        % function evaluation step
    err = toll + 1;                 % error
    k = 0;                          % iterator
    errVect = [];                   % error vector
    xkVect = [];                    % solutions vector
    xGuessCell = num2cell(xGuess);  % inital guess cell array


    %%% Newton method body
    % Check if initial guess is already a solution
    if f(xGuessCell{:,1}) == zeros(length(xGuess),1)
        elapsedTime = toc(timerStart);
        solution = xGuess;
        converge = true;
        
        if nargout == 3
            info = struct;
                info.iteration      = 0;
                info.solutionVector = solution;
                info.errorVector    = f(xGuessCell{:,1});
                info.timeCost       = elapsedTime;
                info.fevalCost      = dimSys;
                info.methodUsed     = 'none';
        end

        % output print on command window
        if visualConfig.print == true
            fprintf('Function zeros are exactly equal to the initial guess\n');
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
            
            switch method
                case 's'
                    % Symbolic math method
                    xkCell = num2cell(xk);          % xk convertion in cell array
                    fk = f(xkCell{:,1});            % function f evaluated in xk
                    feval = feval + dimSys;         % function evaluation counter update
        
                    Jk = J(xkCell{:,1});                % jacobian f evaluated in xk
                    feval = feval + dimSys^(dimSys);    % function evaluation counter update
        
                    methodString = 'SYMBOLIC MATH';  % method used saved for the output

                case 'f'
                    % Finite forward difference method
                    epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));     % epsilon coputed at step xk (vector: dimSys,1)
                    Jk = zeros(dimSys);                             % jacobian initialization
        
                    xkCell = num2cell(xk);                          % xk convertion in cell array
                    fk = f(xkCell{:,1});                            % function f evaluated in xk
                    feval = feval + dimSys;                         % function evaluation counter update
        
                    % loop to evaluate jacobian's i^th column in xk using forward difference
                    for i=1:dimSys
                        epsilonVect = (zeros(dimSys,1));            % epsilon vector initialization
                        epsilonVect(i,1) = epsilon(i,1);            % epsilon vector definition
        
                        xkForward =  xk + epsilonVect;              % xk+epsilon computation in vector
                        xkForwardCell = num2cell(xkForward);        % xk+epsilon convertion in cell array
                        fkForward = f(xkForwardCell{:,1});          % function f evaluated in xk+epsilon
                        feval = feval + dimSys;                     % function evaluation counter update
        
                        Jk(:,i) = (fkForward - fk) / epsilon(i,1);  % i^th columns of jacobian evaluated in xk
                    end
        
                    methodString = 'FORWARD DIFF';                  % method used saved for the output

                case 'c'
                    % Finite centered difference method
                    epsilon = max(sqrt(eps),sqrt(eps)*abs(xk));     % epsilon coputed at step xk (vector: dimSys,1)
                    Jk = zeros(dimSys);                             % jacobian initialization

                    xkCell = num2cell(xk);                          % xk convertion in cell array
                    fk = f(xkCell{:,1});                            % function f evaluated in xk
                    feval = feval + dimSys;                         % function evaluation counter update

                    % loop to evaluate jacobian's i^th column in xk using centered difference
                    for i=1:dimSys
                        epsilonVect = (zeros(dimSys,1));            % epsilon vector initialization
                        epsilonVect(i,1) = epsilon(i,1);            % epsilon vector definition

                        xkForward =  xk + epsilonVect;              % xk+epsilon computation in vector
                        xkForwardCell = num2cell(xkForward);        % xk+epsilon convertion in cell array
                        fkForward = f(xkForwardCell{:,1});          % function f evaluated in xk+epsilon
                        feval = feval + dimSys;                     % function evaluation counter update

                        xkBackward =  xk - epsilonVect;             % xk-epsilon computation in vector
                        xkBackwardCell = num2cell(xkBackward);      % xk-epsilon convertion in cell array
                        fkBackward = f(xkBackwardCell{:,1});        % function f evaluated in xk-epsilon
                        feval = feval + dimSys;                     % function evaluation counter update

                        Jk(:,i) = (fkForward - fkBackward) / (2*epsilon(i,1));   % i^th columns of jacobian evaluated in xk
                    end

                methodString = 'CENTERED DIFF';  % method used saved for the output

            end

            % check if f is sufficiently well-behaved
            if rank(Jk) == zeros(length(xGuess))
                error('Jk has become zero, f not sufficiently well-behaved near solution\n');

            else
                y = Jk \ (-fk);
                xNew = xk + y;              % solutions for the k^th iteration
                err = norm(y);              % error computation at k^th iteration

                errVect = [errVect, err];   % error vector update
                xkVect  = [xkVect, xNew];   % solution vector update
                k = k + 1;                  % iteration index update
                xk = xNew;                  % starting point update for the next iteration
            
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

    if nargout == 3
        info = struct;      % info struct for detailed output
            info.iteration      = k;
            info.solutionVector = xkVect;
            info.errorVector    = errVect;
            info.timeCost       = elapsedTime;
            info.fevalCost      = feval;
            info.methodUsed     = methodString;
    end

    % output print on command window
    if visualConfig.print == true
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
    if visualConfig.plot == true
        iterationVector = (1:k);
        semilogy(iterationVector,errVect,'o-')
    end

end