function [bdfOptions] = bdfSettings(order,startup,method,options,alpha,beta,gamma)
%BDF SETTINGS - Create the options struct for Backward Difference methods
%
%   Syntax:
%       [bdfOptions] = bdfSettings(order,startup,method,options,alpha,beta,gamma)
%
%   Input:
%       order,          double:  select AB method order
%       startup(*),       char:  see startupGuess.m for details
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       options(*),     fsolve:  optimization options, for details see optimoptions.m
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,1]:  insert custom parameter beta for AB method
%       gamma(*),  double[1,n]:  insert custom parameter gamma for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AB method:
%           - amOptions.order = order
%           - amOptions.startup = startup;
%           - amOptions.method = method
%           - amOptions.options = options
%           - amOptions.alpha = alpha
%           - amOptions.beta = beta
%           - amOptions.gamma = gamma
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside ab# functions
%       options: set to don't display iteration by default inside am# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%       gamma:   set as empty by default
%


    %%% Default value for optional input
    if nargin < 7
        gamma = [];
        if nargin < 6    
            beta = [];
            if nargin < 5
                alpha = [];
                if nargin < 4
                    options = [];
                    if nargin < 3
                        method = [];
                        if nargin < 2
                            startup = [];
                        end
                    end
                end
            end
        end
    end


    %%% Struct definition
    bdfOptions = struct;
        bdfOptions.order   = order;
        bdfOptions.startup = startup;
        bdfOptions.method  = method;
        bdfOptions.options = options;
        bdfOptions.alpha   = alpha;
        bdfOptions.beta    = beta;
        bdfOptions.gamma   = gamma;

end