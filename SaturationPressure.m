%Psat1
%Mickey Huang
%Find Vapor Pressure of C3H6 at different temperatures
%Based on Antoine fit of tabulated data
 
 
 
function [VaporPressure] = Psat1(T)
    format long
    %Critical Temperature
    CritTemp = 365.6;
    
    %CASE 1: INVALID T
    %Throw errors and complain a lot if you get a bad (high/low) T input
    if(T<170 || T>CritTemp)
        disp('ERROR: INVALID TEMPERATURE');
        VaporPressure = NaN;
        
    
    %CASE 2: From 170K to BP
    elseif(T>=170 && T<=225.46)
        Y=9.1667-1838.8271/(T-24.5651);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure);
        
    %CASE 3: BP to CP
    elseif(T>225.46 && T<=356.6)
        Y=11.4006-2972.6296/(T+35.8244);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure)
    end
end



%Psat2
%Mickey Huang
%Find Vapor Pressure of 1-C4H8 at different temperatures
%Based on Antoine fit of tabulated data
 
 
 
function [VaporPressure] = Psat2(T)
    format long
    %Critical Temperature
    CritTemp = 418.6;
    
    %CASE 1: INVALID T
    %Throw errors and complain a lot if you get a bad (high/low) T input
    if(T<170 || T>CritTemp)
        disp('ERROR: INVALID TEMPERATURE');
        VaporPressure = NaN;
        
    
    %CASE 2: From 170K to BP
    elseif(T>=170 && T<=266.65)
        Y=9.745-2336.5203/(T-26.614);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure);
        
    %CASE 3: BP to CP
    elseif(T>266.65 && T<=CritTemp)
        Y=9.6916-2407.7014/(T-18.0434);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure)
    end
end

%Psat3
%Mickey Huang
%Find Vapor Pressure of n-C4H10 at different temperatures
%Based on Antoine fit of tabulated data
 
 
 
function [VaporPressure] = Psat3(T)
    format long
    %Critical Temperature
    CritTemp = 425.2;
    
    %CASE 1: INVALID T
    %Throw errors and complain a lot if you get a bad (high/low) T input
    if(T<170 || T>CritTemp)
        disp('ERROR: INVALID TEMPERATURE');
        VaporPressure = NaN;
        
    
    %CASE 2: From 170K to BP
    elseif(T>=170 && T<=272.65)
        Y=9.3851-2292.2293/(T-27.9602);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure);
        
    %CASE 3: BP to CP
    elseif(T>272.65 && T<=CritTemp)
        Y=10.2311-2861.4834/(T+8.2712);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure)
    end
end


%Psat4
%Mickey Huang
%Find Vapor Pressure of cis-2-C4H8 at different temperatures
%Based on Antoine fit of tabulated data
 
 
 
function [VaporPressure] = Psat4(T)
    format long
    %Critical Temperature
    CritTemp = 435.5;
    
    %CASE 1: INVALID T
    %Throw errors and complain a lot if you get a bad (high/low) T input
    if(T<170 || T>CritTemp)
        disp('ERROR: INVALID TEMPERATURE');
        VaporPressure = NaN;
        
    
    %CASE 2: From 170K to BP
    elseif(T>=170 && T<=276.8)
        Y=9.4383-2310.2283/(T-31.5591);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure);
        
    %CASE 3: BP to CP
    elseif(T>272.65 && T<=CritTemp)
        Y=10.0635-2711.8215/(T-6.9010);
        VaporPressure=exp(Y);
        fprintf('%f bars.\n',VaporPressure)
    end
end


