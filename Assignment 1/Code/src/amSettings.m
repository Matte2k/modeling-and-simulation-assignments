function [amOptions] = amSettings(order,startup,method,options,alpha,beta)
%AM SETTINGS - Create the options struct for Adams Moulton methods
%
%   Syntax:
%       [amOptions] = amSettings(order,startup,method,options,alpha,beta)
%
%   Input:
%       order,          double:  select AB method order
%       startup(*),       char:  see startupGuess.m for all the possibility
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       options(*),     fsolve:  optimization options, for details see optimoptions.m
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,n]:  insert custom parameter beta for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AM method:
%           - amOptions.order = order
%           - amOptions.startup = startup;
%           - amOptions.method = method
%           - amOptions.options = options
%           - amOptions.alpha = alpha
%           - amOptions.beta = beta
%
%   Default settings for optional input (*):
%       startup: set as empty by default
%       method:  set as 'Standard' by default inside am# functions
%       options: set to don't display iteration by default inside am# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%


    %%% Default value for optional input
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


    %%% Struct definition
    amOptions = struct;
        amOptions.order   = order;
        amOptions.startup = startup;
        amOptions.method  = method;
        amOptions.options = options;
        amOptions.alpha   = alpha;
        amOptions.beta    = beta;

end
