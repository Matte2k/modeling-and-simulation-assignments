function [tOptions] = tSettings(theta,method,mode,toll,nmax)
%T SETTINGS - Create the options struct for theta method
%
%   Syntax:
%       [tOptions] = tSettings(theta,method,mode,toll,nmax)
%
%   Input:
%       theta(*),       double:  define the 'theta' parameter of the theta methods
%       method(*),        char:  define the solver (available 'fsolve', 'fzero' and 'newton')
%       mode:             char:  set how to compute jacobian (needed only with 'newton')
%       toll:           double:  set tollerance in newton method (needed only with 'newton')
%       nmax:           double:  set max number of iterations in newton method (needed only with 'newton')
%   
%   Output:
%       tOptions,  struct:  contains settings for IEX4 method:
%           - tOptions.theta = theta
%           - tOptions.method = method
%           - tOptions.submethod.mode = mode (present only if method='newton')
%           - tOptions.submethod.toll = toll (present only if method='newton')
%           - tOptions.submethod.nmax = nmax (present only if method='newton')
%
%   Default settings for optional input (*):
%       theta:   set as '0.5' by default inside theta function
%       method:  set as 'fsolve' by default inside theta function
%       mode:    set as default value inside netwon.m
%       toll:    set as default value inside netwon.m
%       nmax:    set as default value inside netwon.m
%


    %%% Default value for optional input
    if nargin < 5
        nmax = [];
        if nargin < 4
            toll = [];
            if nargin < 3
                mode = [];
                if nargin < 2
                    method = [];
                    if nargin < 1
                        theta = [];
                    end
                end
            end
        end
    end


    %%% Struct definition
    tOptions = struct;
        tOptions.theta = theta;
        tOptions.method = method;

        if not(isequal(method,'newton')) 
            tOptions.submethod = [];
            
            if not(isempty(mode)) || not(isempty(toll)) || not(isempty(nmax))
                warning('some input may be unused for the current method selected')
            end
        else
            tOptions.submethod = struct;
                tOptions.submethod.mode = mode;
                tOptions.submethod.toll = toll;
                tOptions.submethod.nmax = nmax;
        end
        
end