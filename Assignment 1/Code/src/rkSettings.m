function [rkOptions] = rkSettings(order,method,alpha,beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rkOptions = struct;
    rkOptions.order = order;
    rkOptions.method = method;
    rkOptions.alpha = alpha;
    rkOptions.beta = beta;
end

