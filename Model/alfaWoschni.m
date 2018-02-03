function [alfa]=alfaWoschni(reducedCa,T,p,pm,Tr,pr,Vr)
global SOC EOC Stroke omega VDisp Bore

C1 = 6.18;                              % Intake / Exhaust period
C2 = 0;

if reducedCa > -180 && reducedCa < SOC  % Compression period
 C1 = 2.28;
 C2 = 0;
end
if reducedCa > SOC && reducedCa < EOC   % Combustion period
 C1 = 2.28;
 C2 = 3.24e-3;
end

Sp  = 2 * Stroke * omega / 2 / pi;
fac = VDisp*Tr/pr/Vr;     

w   = C1*Sp+C2*fac*(p-pm);
if (w<0)
   w=C1*Sp; 
end

C = 3.26;
m = 0.8;
alfa = C*Bore^(m-1)*(p/1000)^m*w^m*T^(0.75-1.62*m);

