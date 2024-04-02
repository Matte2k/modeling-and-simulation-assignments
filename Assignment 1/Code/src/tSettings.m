function [tOptions] = tSettings(theta,method,submethod)
%TSETTINGS Summary of this function goes here
%   Detailed explanation goes here

if isempty(submethod)
    submethod.mode = [];
    submethod.toll = [];
    submethod.nmax = [];
end

tOptions = struct;
    tOptions.theta = theta;
    tOptions.method = method;
    tOptions.submethod = struct;
        tOptions.submethod.mode = submethod.mode;
        tOptions.submethod.toll = submethod.toll;
        tOptions.submethod.nmax = submethod.nmax;
        
end

