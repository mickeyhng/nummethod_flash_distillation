function [x y V L P] = MHflash(zF, T ,P, L, y1)
%% Input check    
    %zF check
    m=iscolumn(zF);
    if m ~= 1
        disp('zF must be a column vector');
        x=NaN;
        y=NaN;
        V=NaN;
        L=NaN;
        P=NaN;
        return
    end
    n=size(zF);
    if n ~= 4
        disp('4 input stream required');
        x=NaN;
        y=NaN;
        V=NaN;
        L=NaN;
        P=NaN;
        return
    end
    F=sum(zF);
    z=zF/F;
    %To make sure the mole fraction add up to 1 to prevent negative stream
    if sum(z)~=1
        disp('invalid stream')
        x=NaN;
        y=NaN;
        V=NaN;
        L=NaN;
        P=NaN;
        return
    end
    %Temperature check
    Psat = [Psat1(T) Psat2(T) Psat3(T) Psat4(T)]';
    if(sum(isnan(Psat))>0)
        disp('ERROR: INVALID TEMPERATURE');
        x=NaN;
        y=NaN;
        V=NaN;
        L=NaN;
        P=NaN;
        return
    end
    if size(T)~=1
        disp('invalid temperature dimension')
        x=NaN;
        y=NaN;
        V=NaN;
        L=NaN;
        P=NaN;
        return
    end
    
    %Valid P --> run flash1.
    if(P>0)
        %Pass arguments to flash1 subroutine.
        [x y V L] = flash1(zF, T, P);
        if(sum(isnan(cat(1,x,y,V,L)))<1) %if any input invalid, isnan(invalid input)=1, error
            %%Pacakage into an array cat function % this step check error
            return
        end
    end
    %valid L --> run flash2.
    if(L>0 && L<sum(zF)) %If L is in the domain of (0, sum(zF))
        %Pass arguments to flash2 subroutine.
        [x y V P] = flash2(zF, T, L);
        if(sum(isnan(cat(1,x,y,V,P)))<1) %exit program if no error
            return
        end
    end
    %valid y1-->run flash3.
    if(y1>0 && y1<1)
        %Pass arguments to flash3 subroutine.
        [x y2 y3 y4 V L P] = flash3(zF, T, y1);
        if(sum(isnan(cat(1,x,y2,y3,y4,V,L,P)))<1)
            %exit program if no error
            return
        end
    %If all three flash fails
    else
        disp('Your conditions are beyond the range of this program, please update');
        x=NaN;
        y=NaN;
        V=NaN;
        L=NaN;
        P=NaN;
        return
    end
end
