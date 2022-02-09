function [x y V P] = flash2(zF, T, L)
disp('Temperature in Kelvin');
disp('zF must be entered in the order of Stream B, E, F, G to get the right answer');
%% input check   
%set maximum number of iterations for the Newton's method algorithm
    kmax = 10000;
    epsilon = 1e-6;
    %Calculate Vapor Pressures at the input Temperature and pressure
    m=iscolumn(zF);
        
    if m ~= 1
        disp('zF must be a column vector');
        x = NaN;
        y = NaN;
        V = NaN;
        P = NaN;
        return
    end
    n=size(zF);
    if n ~= 4
        disp('4 input stream required');
        x = NaN;
        y = NaN;
        V = NaN;
        P = NaN;
        return
    end
    %Calculate Vapor Pressures at the input Temperature
    Psat = [Psat1(T) Psat2(T) Psat3(T) Psat4(T)]';
    %Check that the input temperature is able to return a vapor pressure,
    %otherwise throw an error.
    if(isnan(Psat)==1)
        x = NaN;
        y = NaN;
        V = NaN;
        P = NaN;
        return
    end
    if  size(L)~=1
        disp('L dimension incorrect')
        x = NaN;
        y = NaN;
        V = NaN;
        P = NaN;
        return
    end
    if L <= 0
        disp('No L detected');
        x = NaN;
        y = NaN;
        V = NaN;
        P = NaN;
        return
    end
    if size(T) ~= 1
        disp('T dimension incorrect')
        x = NaN;
        y = NaN;
        V = NaN;
        P = NaN;
        return
    end
    F = sum(zF);
    z = zF/F;
    if sum(z) ~= 1
        disp('mole fraction does not add up')
        return
    end
    %% Initialization
    
    %Initial Guesses for Newton's method
    x{1} = z;
    y{1} = z;
    V{1} = F-L;
    P{1} = 1;
    
    %all variables packaged into an array with vertcat function
    A{1} = vertcat(x{1},y{1},V{1},P{1});
    %% Newton Raphson
    %Begin iteration
    for i=1:kmax
        %Jacobian Matrix
        %x1...x4 y1...y4 V P
           
        J = [1 1 1 1 0 0 0 0 0 0; %sumx
            0 0 0 0 1 1 1 1 0 0; %sum y
            -Psat(1) 0 0 0 P{i} 0 0 0 0 y{i}(1); %Py=Psatx
            0 -Psat(2) 0 0 0 P{i} 0 0 0 y{i}(2);
            0 0 -Psat(3) 0 0 0 P{i} 0 0 y{i}(3);
            0 0 0 -Psat(4) 0 0 0 P{i} 0 y{i}(4);
            L 0 0 0 V{i} 0 0 0 y{i}(1) 0; %zF=yV+xL
            0 L 0 0 0 V{i} 0 0 y{i}(2) 0;
            0 0 L 0 0 0 V{i} 0 y{i}(3) 0;
            0 0 0 L 0 0 0 V{i} y{i}(4) 0];
        %residuals
        r(1) = sum(x{i}) - 1;
        r(2) = sum(y{i}) - 1;
        r(3:6) = P{i}*y{i}-Psat.*x{i};
        r(7:10) = y{i}*V{i}+x{i}*L-zF;
        del=-J\r';
        %Take a newton step
        A{i+1} = A{i}+del;
        %repackage variables into vars array
        x{i+1} = A{i+1}(1:4);
        y{i+1} = A{i+1}(5:8);
        V{i+1} = A{i+1}(9);
        P{i+1} = A{i+1}(10);
        %convergence
       if(abs(norm(A{i+1})-norm(A{i}))/norm(A{i})<epsilon)
            x=x{end}
            y=y{end}
            V=V{end}
            P=P{end}
            %Check condition
            if(V/F<=0 || V/F>=1)
                x = NaN;
                y = NaN;
                V = NaN;
                P = NaN;
                disp('tempearture-Pressure condtion outside of bubble-dew range');
                return
            end
            return;
        end
    end
    disp('Max iteration reached, did not reach convergence')
                x = NaN;
                y = NaN;
                V = NaN;
                P = NaN;
        return
end
