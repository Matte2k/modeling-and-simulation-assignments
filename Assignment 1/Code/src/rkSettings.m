function [rkOptions] = rkSettings(order,method,alpha,beta,iterations)
%RK SETTINGS - Create the options struct for Runge-Kutta methods
%
%   Syntax:
%       [rkOptions] = rkSettings(order,method,alpha,beta)
%
%   Input:
%       order,          double:  select RK method order
%       method(*),        char:  chose between 'Standard' parameters and 'Custom'
%       alpha(*),  double[1,1]:  insert custom parameter alpha for RK method
%       beta(*),   double[1,n]:  insert custom parameter beta for RK method
%       iterations(*),  double:  number of iterations to compute average CPU time
%
%   Output:
%       rkOptions,  struct:  contains settings for RK method:
%           - rkOptions.order = order
%           - rkOptions.method = method
%           - rkOptions.alpha = alpha
%           - rkOptions.beta = beta
%           - rkOptions.iterations = iterations
%
%   Default settings for optional input (*):
%       method:  set as 'Standard' by default inside rk# functions
%       alpha:   set as empty by default
%       beta:    set as empty by default
%       iterations: set as 1 by default
%

    %%% Default value for optional input
    if nargin < 5
        iterations = 1;
        if nargin < 4
            beta = [];
            if nargin < 3
                alpha = [];
                if nargin < 2
                    method = [];
                end
            end
        end
    end


    %%% Struct definition
    rkOptions = struct;
        rkOptions.order  = order;
        rkOptions.method = method;
        rkOptions.alpha  = alpha;
        rkOptions.beta   = beta;
        rkOptions.iterations = iterations;

end

