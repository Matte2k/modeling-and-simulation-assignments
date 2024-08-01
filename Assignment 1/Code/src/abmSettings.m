function [abmOptions] = abmSettings(order,startup,method,alpha,betaP,betaC)
%AB MSETTINGS - Create the options struct for Adams Bashforth Moulton methods
%
%   Syntax:
%       [abmOptions] = abmSettings(order,startup,method,alpha,betaP,betaC)
%
%   Input:
%       order,          double:  select ABM method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom' 
%       alpha(*),  double[1,1]:  insert custom parameter alpha for ABM method
%       betaP(*),  double[1,n]:  insert custom parameter beta predictor for ABM method
%       betaP(*),  double[1,n]:  insert custom parameter beta corrector for ABM method
%
%   Output:
%       abmOptions,  struct:  contains settings for AB method:
%           - abmOptions.order = order
%           - abmOptions.startup = startup;
%           - abmOptions.method = method
%           - abmOptions.alpha = alpha
%           - abmOptions.betaP = betaP
%           - abmOptions.betaC = betaC
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside abm# function
%       alpha:   set as empty by default
%       betaP:    set as empty by default
%       betaC:    set as empty by default
%


    %%% Default value for optional input
    if nargin < 6    
        betaC = [];
        if nargin < 5
            betaP = [];
            if nargin < 4          
                alpha = [];
                if nargin < 3
                    method = [];
                    if nargin < 2
                        startup = [];
                    end
                end
            end
        end
    end

    %%% Struct definition
    abmOptions = struct;
        abmOptions.order   = order;
        abmOptions.startup = startup;
        abmOptions.method  = method;
        abmOptions.alpha   = alpha;
        abmOptions.betaP   = betaP;
        abmOptions.betaC   = betaC;

end