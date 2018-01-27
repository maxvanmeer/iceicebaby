function [alfa]=alfaWoschni(reducedCa,T,p,pm,Vd,Tr,pr,Vr)
global SOC
C1 = 6.18;
C2 = 0;
if reducedCA > -180 && reducedCa < SOC % Compression period
 C1 = 2.28;
 C2 = 0;
end
if reducedCA > SOC && reducedCa < 180 % Compression period
 C1 = 2.28;
 C2 = 3.24e-3;
end
fac = Vd*Tr/pr/Vr;
w = C1*Sp+C2*fac*(p-pm);
alfa = C*B^(m-1)*p^m*w^m*T^(0.75-1.62*m);