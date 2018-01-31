function [yp] = FtyDAE( t,y )
global Int Exh QLHV SpS Runiv omega
global GaussatCA50  mfuIVCClose si EtaComb Bore Stroke Omega rc CA50 BDUR iCase

Twall   = 273+80;   %80 degrees Celcius is on the lower side. Higher is better for the engine efficieny. 
Tpiston = 273+110;  %110 degrees Celsius is a guess, based on the normal temperature of engine oil (which cools the pistons). NOT SURE.
alfa    = 500;

VDisp   = pi*(Bore/2)^2*Stroke;
Vc      = VDisp/(rc-1);

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Mi = [SpS.Mass];Nsp = length(Mi);
%%
pI=Int.p;TI=Int.T;
pE=Exh.p;TE=Exh.T;
%%
p=y(1);T=y(2);mi=y(3:end); %Note that p and T are swapped compared to the slides
m = sum(mi);
[V,dVdt,A,A_c,A_p]=CylVolumeFie(t);
Yi = [mi/m]';
%%
Ca          = time2ang(t,omega)/(2*pi)*360;
reducedCa   = mod(Ca+360,720)-360;
CADS        = omega/(2*pi)*360;
%%
for ii=1:Nsp
    hi(ii) = HNasa(T,SpS(ii));
    ei(ii) = ENasa(T,SpS(ii));
    Cpi(ii) = CpNasa(T,SpS(ii));
    Cvi(ii) = CvNasa(T,SpS(ii));
end
h = hi*Yi';
e = ei*Yi';
Cp = Cpi*Yi';
Cv = Cvi*Yi';
gamma = Cp/Cv;
%% Intake and Exhaust thermodynamical properties
hI      = [Int.h];
eI      = [Int.e];
CpiI    = [Int.Cp];
hE      = [Exh.h];
eE      = [Exh.e];
CpiE    = [Exh.Cp];
%%
StateCyl.p          = p;
StateCyl.T          = T;
StateCyl.gamma      = gamma;
StateCyl.rho        = m/V;
StateCyl.Rg         = Cp-Cv;

StateIntake.p       = pI;
StateIntake.T       = TI;
YY                  = [Int.Y];
Ma                  = 1/(sum(YY./Mi));
Rg                  = Runiv/Ma;
StateIntake.rho     = pI/TI/Rg;
hpI                 = YY*hI';
CpI                 = YY*CpiI';
StateIntake.gamma   = CpI/(CpI-Rg);
[dmdtI]             = mdot(StateIntake,StateCyl,t,Int);


StateExhaust.p      = pE;
StateExhaust.T      = TE;
YY                  = [Exh.Y];
Ma                  = 1/(sum(YY./Mi));
Rg                  = Runiv/Ma;
hpE                 = YY*hE';       
StateExhaust.rho    = pE/TE/Rg;
CpE                 = YY*CpiE';
StateExhaust.gamma   = CpE/(CpE-Rg);
[dmdtE]             = mdot(StateExhaust,StateCyl,t,Exh);

if (abs(dmdtI) > 0)
    mfuIVCClose = mi(1);
end

if (dmdtI > 0)
    hpaI = hpI;
    YI = Int.Y;
else
    hpaI = h;
    YI = Yi;
end
if (dmdtE > 0)
    YE    = Yi;     % Same as cylinder but enthalpy is that of exhaust
    hpaE  = YE*hE';
else
    hpaE  = h;
    YE    = Yi;
    Exh.Y = YE;
end
dmidt       = [YI*dmdtI + YE*dmdtE]';
%dmfuComb    = EtaComb*mfuIVCClose*pdf(GaussatCA50,reducedCa)*CADS;
dmfuComb    = EtaComb*mfuIVCClose*wiebefunctions(reducedCa)*CADS;
dmidt_c     = si'*dmfuComb;
dmidt       = dmidt - dmidt_c;
dQcomb      = QLHV*dmfuComb;
dQcomb_real = ei*dmidt_c;

%% Woschni heatloss
%normal model
% CA50=10;                        % CA50 (50% HR), replace with 
load('currentCase.mat');
HR = currentCase.HR;
ReducedCA = -360:360;

CA10 = ReducedCA(find(HR>0.1*HR(length(HR)),1)+1);%Ze liggen een achter omdat ReducedCA van 0 tot 360 loopt en HR van 1 tot 360
CA50 = ReducedCA(find(HR>0.5*HR(length(HR)),1)+1);
CA90 = ReducedCA(find(HR>0.9*HR(length(HR)),1)+1);
% BDUR=20;                        % Burn Duration, replace with data case
BDUR = CA90-CA10;
CAign = CA50-0.5*BDUR;
CAend = CAign+BDUR;
SOC = CAign;
EOC = CAend;

V0      = CylVolumeFie(t(1));
T0      = 273;
p0      = 3.5*10^5;                                             %might need to be adjusted because of looad change

pm = ((VDisp+Vc)/(V))^gamma * p0;      

alfa = alfaWoschni(reducedCa,T,p,pm,T0,p0,V0);

dQhl = alfa*(A_c*(Twall-T)+A_p*(Tpiston-T));

%% DAE formulation
Rg = StateCyl.Rg;
yp = [dQhl-p*dVdt+hpaI*dmdtI+hpaE*dmdtE;
    p*V-Rg*T*m;
    dmidt];
end

