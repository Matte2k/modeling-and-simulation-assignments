function [F,guess,degVec,hVec] = selectOp(method,dimSys,params)
%SELECT OP - Finds the F(h,alpha) operator and associated settings usefull
% to compute stability region of the desired integration method 
%
%    Syntax:
%       [F,guess,degVec,hVec] = selectOp(method,dimSys,params)
%
%   Input:
%       method,       char:  select method to obtain the operator
%       dimSys(*),  double:  analyzed system dimension 
%       params(*),  double:  input used in BI2 method (sets theta value)
%
%   Output:
%       F,      function(h,A):  corresponding method operator
%       guess,         double:  eductaed initial guess used for fzero in stabRegion.m
%       degVec,   double[1,2]:  starting and ending alpha used in stability region analysis
%       hVec,     double[1,2]:  starting and ending h used in stability guess analysis
%
%   Default settings for optional input (*):
%       dimSys:     set to 2 by default for stability region analysis
%       params:     set to 0.4 by default
%
%   Note:   'guess','degVec','hVec' has been pre defined throught an iterarive
%           based on graphic analysis obtained using the function stabGuess.m
%


%%% Optional input definition
if nargin < 2
    dimSys = 2;
    if nargin < 3       
        params = 0.4;
    end
end
         

%%% Operator selection
I = eye(dimSys);
switch method
    case 'RK1'
        F = @(h,A) I + (h*A);
        guess = 2;
        degVec = [180 0];
        hVec = [0 5];

    case 'RK2'
        F = @(h,A) I + (h*A) + 1/2*(h*A)^2;
        guess = 2;
        degVec = [180 0];
        hVec = [0 5];

    case 'RK3'
        F = @(h,A) I + (h*A) + 1/2*(h*A)^2 + 1/6*(h*A)^3;
        guess = 2;
        degVec = [180 0];
        hVec = [0 5];

    case 'RK4'
        F = @(h,A) I + (h*A) + 0.5*(h*A)^2 + 1/6*(h*A)^3 + 1/24*(h*A)^4;
        guess = 2;
        degVec = [180 0];
        hVec = [0 5];

    case 'BI2'
        F = @(h,A) (I-(1-params)*h*A+((1-params)*h*A)^2/2) \ (I+params*h*A+(params*h*A)^2/2);

        if params < 0.5
            guess  = 10;
            degVec = [0 180];
            hVec = [0 12];

        elseif params >= 0.5    % DEBUG params = 0.5
            guess  = 10;
            degVec = [180 0];
            hVec = [0 12];
        end
        
    case 'IEX4'
        F  = @(h,A) -1/6*((I-h*A/1)\I)^1 + 4*((I-h*A/2)\I)^2 - 27/2*((I-h*A/3)\I)^3 + 32/3*((I-h*A/4)\I)^4;
        guess  = 10;
        degVec = [0 180];
        hVec = [10 15];
           
%    % TO BE FIXED the next ones
%    case 'AB1'
%        %F = @(h,A) (h*A + I) / (1);
%    case 'AB2'
%        %F = @(h,A) 2*(h*A^2-h*A)/(3*h*A-1);
%        %F = @(h,A) ((h*A)^2 - (h*A)) / ((3*(h*A) - 1)/2);
%        %F = @(h,A) (I + (h*A))^1 + h/2 * ( 3 * ((I + (h*A))^1) * A - A );
%    case 'AB3'
%        %F = @(h,A) ((h*A)^3 - (h*A)^2) / ((5 - 16*(h*A) + 23*(h*A)^2)/12);
%    case 'AB4'
%        %F = @(h,A) ((h*A)^4 - (h*A)^3) / ((55*(h*A)^3 - 59*(h*A)^2 + 37*(h*A) - 9)/24);
%    case 'AM1'
%    case 'AM2'
%    case 'AM3'
%    case 'AM4'
%    case 'ABM1'
%    case 'ABM2'
%    case 'ABM3'
%    case 'ABM4'
%    case 'BDF1'
%        %F = @(h,A) (1 - 1/(h*A));
%    case 'BDF2'
%        %F = @(h,A) (1 - 1/(h*A)) + (1 - 1/(h*A))^2 / 2;
%    case 'BDF3'
%        %F = @(h,A) (1 - 1/(h*A)) + (1 - 1/(h*A))^2 / 2 + (1 - 1/(h*A))^3 / 3;
%    case 'BDF4'
%        %F = @(h,A) (1 - 1/(h*A)) + (1 - 1/(h*A))^2 / 2 + (1 - 1/(h*A))^3 / 3 + (1 - 1/(h*A))^4 / 4;

    otherwise
        error('Insert a valid method')

end

end