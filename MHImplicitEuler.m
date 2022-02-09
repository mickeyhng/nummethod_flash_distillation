%Implicit Euler method to solve an unsteady state diffusion problem with no
%chemical reaction
%The algorithm used here is derived from the theory lay out 
%https://opencommons.uconn.edu/cgi/viewcontent.cgi?referer=https%3A%2F%2Fwww.google.com%2F&httpsredir=1&article=1118&context=srhonors_theses
%The program will run 1500 iterations
function MHImplicitEuler
%% Initialization
D=3.4e-5; % m^2/2 diffusivity of methane in air at 400K 
di=0; %initial distance m
df=1; %final distance m
ti=0; %Time initial from 0 to 10 seconds specifed in the problem
tf=10; 
hd=(df-di)/1500;% distance step size
ht=(tf-ti)/1500; %Stepsize for time
x=linspace(di,df,1501);
t=linspace(ti,tf,1501); %creating equally spacing vector
dL=hd; %delx and delt, the smaller the step size, the more accurate the solution is but at a price of higher computational cost
dt=ht;
%% setting up the and backward differencing
% a banded matrix was set up in the source guide on page 9, this matrix is
% used to solve for the implicit solution of this parabolic PDE
a=D*dt/(dL^2); 
A=zeros(1500,1500); %preallocating space
for i = 1:1500
    A(i,i) = 1 + 2*a; %diagonal matrix
end
for i = 1:1499     
    A(i,i+1)=-a; %banded matrix
    A(i+1,i)=-a;
end
%% BC dC/dx@x=L = 0 C(0,t)=18 C(x,0)=0
BC = zeros(1,1500)';%preallocating space
BC(1,1)=18; %concentration of methane if there is a leak
F=zeros(1500, 1501);
for i= 1:1500
    F(:,i+1)= A\(BC + F(:,i));
end
for i = 1:1500
    if F(i,end) <= 0.8 
        d = i;    %F is keep getting added for 1500 times, so the L the program will give is the sum distance calculated by all iterations, to find the average, divide the sum by 1500
        b=d/1500; %divide by number of iteration
        fprintf ('%.4f\n', b);
        break
    end
end
