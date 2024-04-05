function [rkOptions] = rkSettings(order,method,alpha,beta)
%RK SETTINGS - Create the options struct for Runge-Kutta methods
%
%   Syntax:
%       [rkOptions] = rkSettings(order,method,alpha,beta)
%
%   Input:
%       order,          double:  select AB method order
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       alpha(*),  double[1,1]:  insert custom parameter alpha for AB method
%       beta(*),   double[1,n]:  insert custom parameter beta for AB method
%
%   Output:
%       abOptions,  struct:  contains settings for AM method:
%           - amOptions.order = order
%           - amOptions.method = method
%           - amOptions.alpha = alpha
%           - amOptions.beta = beta
%
%   Default settings for optional input (*):
%       method:  set as 'Standard' by default inside ab# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%

    %%% Default value for optional input
    if nargin < 4
        beta = [];
        if nargin < 3
            alpha = [];
            if nargin < 2
                method = [];
            end
        end
    end


    %%% Struct definition
    rkOptions = struct;
        rkOptions.order  = order;
        rkOptions.method = method;
        rkOptions.alpha  = alpha;
        rkOptions.beta   = beta;

end

