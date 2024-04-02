function [iOptions] = iSettings(method,options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    if nargin < 1
        method = [];
    end
    options = [];
end

iOptions = struct;
    iOptions.method = method;
    iOptions.options = options;
end

