function [abOptions] = abSettings(order,startup,method,alpha,beta)
%AB SETTINGS - Create the options struct for Adams Bashforth methods
%
%   Syntax:
%       [abOptions] = abSettings(order,startup,method,alpha,beta)
%
%   Input:
%       order,          double:  select AB method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom' 
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,n]:  insert custom parameter beta for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AB method:
%           - abOptions.order = order
%           - abOptions.startup = startup;
%           - abOptions.method = method
%           - abOptions.alpha = alpha
%           - abOptions.beta = beta
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside ab# function
%       alpha:   set as empty by default
%       beta:    set as empty by default
%


    %%% Default value for optional input
    if nargin < 5
        beta = [];
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

    %%% Struct definition
    abOptions = struct;
        abOptions.order   = order;
        abOptions.method  = method;
        abOptions.alpha   = alpha;
        abOptions.beta    = beta;
        abOptions.startup = startup;

end