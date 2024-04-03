function [bdfOptions] = bdfSettings(order,method,options,alpha,beta,gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

bdfOptions = struct;
    bdfOptions.order = order;
    bdfOptions.method = method;
    bdfOptions.options = options;
    bdfOptions.alpha = alpha;
    bdfOptions.beta = beta;
    bdfOptions.gamma = gamma;
end

