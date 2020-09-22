function [x,fval]= unimin(f,deltaM,ver,interval)
if(~exist('deltaM','var'))
    deltaM = 1e-8;
end

if(~exist('ver','var'))
    ver = 0;
end

if(~exist('interval','var'))
    interval = [0 1];
end
phi = (1+sqrt(5))/2;
%phi=2;
x(1) = interval(1);
x(4) = interval(2);
try
fval(1) = f(x(1));
fval(4) = f(x(4));
catch e
   
keyboard 
end

x(2) = x(4) - (x(4)-x(1))/phi;
x(3) = x(1) + (x(4)-x(1))/phi;
fval(2) = f(x(2));
fval(3) = f(x(3));
k=0;
while abs(x(1)-x(4))>deltaM
    k=k+1;
    if(ver)
        fprintf('x:\t %5.3f %5.3f %5.3f %5.3f\n',x)
        fprintf('fx:\t %5.3f %5.3f %5.3f %5.3f\n\n',fval)
    end
    if(fval(2) <fval(3))
        x(4) = x(3);
        fval(4) = fval(3);
        
        x(3) = x(2);
        fval(3) = fval(2);
        
        x(2) = x(4) - (x(4)-x(1))/phi;
        fval(2) =f(x(2));
    else
        x(1) = x(2);
        fval(1) = fval(2);
        
        x(2) = x(3);
        fval(2) = fval(3);
        
        x(3) = x(1) + (x(4)-x(1))/phi;
        fval(3) = f(x(3));
    end
end

[fval,im] = min(fval);
x=  x(im);
end