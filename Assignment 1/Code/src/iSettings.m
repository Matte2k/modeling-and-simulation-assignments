function [iOptions] = iSettings(method,options)
%I SETTINGS - Create the options struct for iex4 methods
%
%   Syntax:
%       [iOptions] = iSettings(method,options)
%
%   Input:
%       method(*),        char:  define the solver (available 'fsolve' only in the current realese)
%       options(*),     fsolve:  optimization options, for details see optimoptions.m
%
%   Output:
%       iOptions,  struct:  contains settings for IEX4 method:
%           - iOptions.method = method
%           - iOptions.options = options
%
%   Default settings for optional input (*):
%       method:  set as 'fsolve' by default inside iex4 function
%       options: set to don't display iteration by default inside iex4 function
%


    %%% Default value for optional input
    if nargin < 2
        options = [];
        if nargin < 1
            method = [];
        end
    end


    %%% Struct definition
    iOptions = struct;
        iOptions.method  = method;
        iOptions.options = options;

end