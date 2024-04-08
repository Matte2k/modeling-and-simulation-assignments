function[F]=selectOp(method,dimSys,params)

% Maybe add starting point suggestion?

if nargin < 3       % Temporary
    params = 0.4;
end
         
I = eye(dimSys);

switch method
    case 'RK1'
        %F = @(h,A) I + (h*A);
                F = @(h,A) I + (h.*A);

        % start in deg=180, guess 2

    case 'RK2'
        %F = @(h,A) I + (h*A) + 1/2*(h*A)^2;
                F = @(h,A) I + (h*A) + 1/2*(h*A)^2;

        % start in deg=180, guess 2

    case 'RK3'
        %F = @(h,A) I + (h*A) + 1/2*(h*A)^2 + 1/6*(h*A)^3;
                F = @(h,A) I + (h*A) + 1/2*(h*A)^2 + 1/6*(h*A)^3;

        % start in deg=180, guess 2

    case 'RK4'
        %F = @(h,A) I + (h*A) + 0.5*(h*A)^2 + 1/6*(h*A)^3 + 1/24*(h*A)^4;
                F = @(h,A) I + (h*A) + 0.5*(h*A)^2 + 1/6*(h*A)^3 + 1/24*(h*A)^4;

        % start in deg=180, guess 2

    case 'BI2'
        theta = params;
        F = @(h,A) (I-(1-theta)*h*A+((1-theta)*h*A)^2/2) \ (I+theta*h*A+(theta*h*A)^2/2);
        % start in deg=0,   guess 10 -> for theta < 0.5
        % start in deg=180, guess 10 -> for theta > 0.5

    case 'IEX4'
        F  = @(h,A) -1/6*((I-h*A/1)\I)^1 + 4*((I-h*A/2)\I)^2 - 27/2*((I-h*A/3)\I)^3 + 32/3*((I-h*A/4)\I)^4;
        % start in deg=0, guess = 10
    
        
        % TO BE FIXED the next ones
    case 'AB1'
        F = @(h,A) (h*A - I) / (1);
    case 'AB2'
        F = @(h,A) (h*A - I) / ((3 - 1/(h*A))/2);
    case 'AB3'
        F = @(h,A) (h*A - I) / ((23 - 16/(h*A) + 5/(h*A)^2) / 12);
    case 'AB4'
        F = @(h,A) (h*A - I) / ((55 - 59/(h*A) + 37/(h*A)^2 - 9/(h*A)^3) / 24);
    case 'AM1'
    case 'AM2'
    case 'AM3'
    case 'AM4'
    case 'ABM1'
    case 'ABM2'
    case 'ABM3'
    case 'ABM4'
    case 'BDF1'
        F = @(h,A) (1 - 1/(h*A));
    case 'BDF2'
        F = @(h,A) (1 - 1/(h*A)) + (1 - 1/(h*A))^2 / 2;
    case 'BDF3'
        F = @(h,A) (1 - 1/(h*A)) + (1 - 1/(h*A))^2 / 2 + (1 - 1/(h*A))^3 / 3;
    case 'BDF4'
        F = @(h,A) (1 - 1/(h*A)) + (1 - 1/(h*A))^2 / 2 + (1 - 1/(h*A))^3 / 3 + (1 - 1/(h*A))^4 / 4;

    otherwise
        error('Insert a valid method')

end


