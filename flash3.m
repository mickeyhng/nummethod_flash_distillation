%%Steepest descent is very senstitive to initial guesses. This algorithm
%%attempts to improve the previous steepest descent algorithm by
%%incorporating Newton's Method to find a better initial guess.
%%This method requires high computational cost, to minimize the cost, the
%%10 equations are reduced to 5 equations and different naming was used
 
function [x y2 y3 y4 V L P] = flash3(zF,T,y1)
Go = 1;
 
%% input check
m=iscolumn(zF);
    if m ~= 1
        disp('zF must be a column vector');
        x = NaN;
        y2 = NaN;
        y3 = NaN;
        y4 = NaN;
        V = NaN;
        P = NaN;
        L = NaN;
        return
    end
    n=size(zF);
    if n ~= 4
        disp('4 input stream required');
        x = NaN;
        y2 = NaN;
        y3 = NaN;
        y4 = NaN;
        V = NaN;
        P = NaN;
        L = NaN;
        return
    end
    %Calculate Vapor Pressures at the input Temperature
    P1=Psat1(T);
    P2=Psat2(T);
    P3=Psat3(T);
    P4=Psat4(T);
    Psat = [P1 P2 P3 P4]';
    
    %Check that the input temperature is able to return a vapor pressure,
    %otherwise throw an error.
    if(isnan(Psat)==1)
        x = NaN;
        y2 = NaN;
        y3 = NaN;
        y4 = NaN;
        V = NaN;
        P = NaN;
        L = NaN;
        return
    end
    if  size(y1)~=1
        disp('y1 dimension incorrect')
        x = NaN;
        y2 = NaN;
        y3 = NaN;
        y4 = NaN;
        V = NaN;
        P = NaN;
        L = NaN;
        return
    end
    if y1 <= 0 || y1>1
        disp('No y1 detected');
        x = NaN;
        y2 = NaN;
        y3 = NaN;
        y4 = NaN;
        V = NaN;
        P = NaN;
        L = NaN;
        return
    end
    if size(T) ~= 1
        disp('T dimension incorrect')
        x = NaN;
        y2 = NaN;
        y3 = NaN;
        y4 = NaN;
        V = NaN;
        P = NaN;
        L = NaN;
        return
    end
    F = sum(zF);
    z = zF/F;
    if sum(z) ~= 1
        disp('mole fraction does not add up')
        return
    end
 
%% Initialization 
% x(1) = L
% x(2) = y2
% x(3) = y3
% x(4) = y4
% x(5) = P
 
f1 = @(x) zF(1)-y1*(F-x(1))-y1*x(5)*x(1)/P1;% f1 = z1*F - y1*(F-L) - y1*P*L/P1
f2 = @(x) zF(2)-x(2)*(F-x(1))-x(2)*x(5)*x(1)/P2;% f2 = z2*F - y2*(F-L) - y2*P*L/P2
f3 = @(x) zF(3)-x(3)*(F-x(1))-x(3)*x(5)*x(1)/P3;% f3 = z3*F - y3*(F-L) - y3*P*L/P3
f4 = @(x) zF(4)-x(4)*(F-x(1))-x(4)*x(5)*x(1)/P4; % f4 = z4*F - y4*(F-L) - y4*P*L/P4
f5 = @(x) y1+x(2)+x(3)+x(4)-1; % f5 = y1 + y2 + y3 + y4 - 1
 
%residual
r = @(x) [f1(x); f2(x); f3(x); f4(x); f5(x)];
%Jacobian
jac = @(x) [y1-y1*x(5)/P1 0 0 0 -y1*x(1)/P1;
           x(2)-x(2)*x(5)/P2 x(5)-x(1)*x(5)/P2 0 0 -x(2)*x(1)/P2;
           x(3)-x(3)*x(5)/P3 0 x(5)-x(1)*x(5)/P3 0 -x(3)*x(1)/P3;
           x(4)-x(4)*x(5)/P4 0 0 x(5)-x(1)*x(5)/P4 -x(4)*x(1)/P4;
           0 1 1 1 0];
%Initial guess
X0 = [F/2 z(2) z(3) z(4) z(1)*P1/y1]';
 
%Newton's method to improve initial guess
   for k = 1:3 % 3 iteratives step recommended by the project guideline. Anything more than that will lead to bad result
    r0 = r(X0);
    jac0 = jac(X0);
    del = -jac0\r0;
    X = X0+del;
    X0 = X;
   end
 
%% Steepest Descent
% Setting up g(x) the square of residual to be used in steepest descent
g = @(x) f1(x)^2+f2(x)^2+f3(x)^2+f4(x)^2+f5(x)^2;
%gradient
delG=@(x) [2*f1(x)*(y1-y1*x(5)/P1)+ 2*f2(x)*(x(2)-x(2)*x(5)/P2)+2*f3(x)*(x(3)-x(3)*x(5)/P3)+2*f4(x)*(x(4)-x(4)*x(5)/P4);
           2*f2(x)*(-(F-x(1))-x(1)*x(5)/P2)+2*f5(x);
           2*f3(x)*(-(F-x(1))-x(1)*x(5)/P3)+2*f5(x);
           2*f4(x)*(-(F-x(1))-x(1)*x(5)/P4)+2*f5(x);
           2*f1(x)*(-y1*x(1)/P1)+2*f2(x)*(-x(2)*x(1)/P2)+2*f3(x)*(-x(3)*x(1)/P3)+2*f4(x)*(-x(4)*x(1)/P4)];
    
 % Descending %the following algorithm is based on armijo's rule on http://www.mit.edu/~dimitrib/PTseng/516/hmwk1_add.txt
for k = 1:20000    
    nf = 0; %nf=number of function evaluation
  for i = 1:100    % for loop to find a
   a = 0.5^nf;   
   d = -delG(X0);
   objfunc = g(X0); 
   newobj = g(X0+a*d);
    if newobj > objfunc % algorithm is not descending, decrease a and check again
       nf = nf + 1;
    else
       Xnext = X0 + a*d; %newobj is descending, break the sub for loop and goes back to main for loop
       break
    end   
  end
  if a*norm(d)/norm(Xnext) < 1e-9 % once convergence, method complete
      break
  else
      X0 = Xnext; %update Xnext to be the next X0
  end
  if i == 100 || k == 20000
      disp('Max iteration reached')  
      Go = 0;
      break
  end   
end
 
% Output value updated
 
L = Xnext(1)
y2 = Xnext(2)
y3 = Xnext(3)
y4 = Xnext(4)
P = Xnext(5)
 
%% Backward calculationn
% Determine x through yi = Kxi = Pi*xi/P
% This rearranges to xi = P*yi/Pi
x=[y1*P/P1;y2*P/P2;y3*P/P3;y4*P/P4]
V=F - L
 
% Checks to make sure that the calculated pressure and temperature
% combinations are valid 
%Very difficult to calculate bubble point temperature
Pbubble = (z(1)*P1+z(2)*P2+z(3)*P3+z(4)*P4);
Pdew = (z(1)/P1+z(2)/P2+z(3)/P3+z(4)/P4)^-1;
if P > Pbubble
    disp('Pressure above bubble-point')
    return
end
if P < Pdew
    disp('Pressure below dew-point')
    return
end
if V/F<=0 || V/F>=1
    disp('Outside of bubble-dew condition')
    return
end
% When Go Marker is turned off
if Go == 0
    x = NaN;
    y2 = NaN;
    y3 = NaN;
    y4 = NaN;
    L = NaN;
    V = NaN;
    P=NaN;
end
