F = ode;
% F.ODEFcn = @(t,y)[ -0.04*y(1) + 10^4 * y(2) * y(3);...
%                     0.04*y(1) - 10^4 * y(2) * y(3) - 3 * 10^7 * y(2)^2;...
%                     3*10^7*y(2)^2 ];
F.ODEFcn = @robertsdae;
F.InitialValue = [1 0 0];
F.AbsoluteTolerance = 1e-12;
F.RelativeTolerance = 1e-6;
F.MassMatrix = [1 0 0; 0 1 0; 0 0 0];
F.Solver = "auto";
sol = solve(F,0,60);
F.SelectedSolver


% y0 = [1; 0; 0];
% tspan = [0 4*logspace(-6,6)];
% M = [1 0 0; 0 1 0; 0 0 0];
% options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',[1e-6 1e-10 1e-6]);
% [t,y] = ode15s(@robertsdae,tspan,y0,options);