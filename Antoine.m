% Use Antoine parameter to determine vapor pressure
%A,B,C Antoine parameter
log(P) = A - B / ( T + C )
log(P) = A - B / ( T + C )
 T.log(P) + C.log(P) = A.T + A.C - B
 log(P) = A + (A.C-B)/T - C.log(P)/T
 %a0 = A a1=A*C-B a2=-C
 %x1=1/T x2=log(P)/T
 y = a0 + a1.x1 + a2.x2
 X = [ones(size(x1)) x1 x2];
 b = regress(y,X) %multiple linear regression return to 3 coefficients
 b(1)=a0=A %A
 -b(3)=a2=C %C
 b(2)=a1
B=A*C-b(2) 

