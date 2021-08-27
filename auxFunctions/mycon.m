function [c,ceq] = mycon(x)
% Enforces constraints for fmincon.
% See fmincon documentation for details.
c = x(3)*2.2 - x(2); 
ceq = 0;