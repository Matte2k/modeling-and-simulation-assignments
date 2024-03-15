function [x,iter,res] = newtonSym(f, x0, toll)
%newtonSym
%   Function that compute the zero of given function using the Newton's
%   method. The jacobian is computed using symbolic manipulation

nSymVars = nargin(f);
fSym = str2sym(char(f));

a=0;

end

