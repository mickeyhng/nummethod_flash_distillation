function [x y V L] = flash1(zF, T, P)
%% input check
disp('Temperature in Kelvin and pressure in bar');
disp('zF must be entered in the order of Stream B, E, F, G to get the right answer');
    %Calculate Vapor Pressures at the input Temperature and pressure
    m=iscolumn(zF);
        
    if m ~= 1
        disp('zF must be a column vector');
        x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
        return
    end
    n=size(zF);
    if n ~= 4
        disp('4 input stream required');
        x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
        return
    end
    if size(P) ~= 1
        disp('P dimension incorrect')
        x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
        return
    end
    if P < 0
        disp('negative input pressure detected');
        x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
        return
    end
    if size(T) ~= 1
        disp('T dimension incorrect')
        x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
        return
    end
    Psat = [Psat1(T) Psat2(T) Psat3(T) Psat4(T)]';
    %Check that the input temperature is able to return a vapor pressure, 
    %otherwise the program will throw an error.
    if(isnan(Psat)==1)
        x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
        return
    end
    %stream properties calculation
    %Approximating equilibrium constant as the ratio of Vapor Pressure to Pressure
    F = sum(zF);
    z = zF/F;
    K = Psat/P;
    %Check that flash is operating within bubble dew point
    Pbubble = (z(1)*Psat1(T)+z(2)*Psat2(T)+z(3)*Psat3(T)+z(4)*Psat4(T));
    Pdew = (z(1)/Psat1(T)+z(2)/Psat2(T)+z(3)/Psat3(T)+z(4)/Psat4(T))^-1;
if P > Pbubble
    disp('Error: Pressure above bubble-point')
    x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
    return
end
if P < Pdew
    disp('Error: Pressure below dew-point')
    x = NaN;
        y = NaN;
        V = NaN;
        L = NaN;
    return
end  
    %Check to see if z adds up to 1, if negative component in zF is
    %detected, then sum(z) can't be equal to 1
    if sum(z)~= 1
        disp('mole fractions of the stream does not add up to 1');
        return
    end
    %% Solving using Rachford-Rice equation
    %Defining Rachford-Rice equation as a function of V, all other variables are known.
    fun= @(V) sum(zF.*(1-K)./(F+V*(K-1)));
    %Solve the equation using root finding algorithm
    V = fzero(fun, sum(zF)); 
    %Calculate other properties of the outlet streams
    L = F - V;
    x = zF./((K-1)*V+F);
    y = K.*x;
    
    L
    V
    xB=x(1)
    xE=x(2)
    xF=x(3)
    xG=x(4)
    yB=y(1)
    yE=y(2)
    yF=y(3)
    yG=y(4)
end
