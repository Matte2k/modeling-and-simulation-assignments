function [abmOptions] = abmSettings(order,method,alpha,betaP,betaC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

abmOptions = struct;
    abmOptions.order = order;
    abmOptions.method = method;
    abmOptions.alpha = alpha;
    abmOptions.betaP = betaP;
    abmOptions.betaC = betaC;

end
