function [yp] = FtyDAE( t,y )
global Int Exh QLHV SpS Runiv omega
global GaussatCA50  mfuIVCClose si EtaComb Bore Stroke Omega rc

Twall   = 273+80;
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
dmfuComb    = EtaComb*mfuIVCClose*pdf(GaussatCA50,reducedCa)*CADS;
dmidt_c     = si'*dmfuComb;
dmidt       = dmidt - dmidt_c;
dQcomb      = QLHV*dmfuComb;
dQcomb_real = ei*dmidt_c;

%% Woshni heatloss
C1 = 2.28;
C2 = 0;
if (dmdtI > 0) || (dmdtE > 0)
C1 = 6.18;
C2 = 0;
phase = "intake/exhaust"
end
if (dmfuComb > 0)
C1 = 2.28;
C2 = 3.24*10-3;
phase = "combustion"
end

V0      = CylVolumeFie(t(1));
T0      = 273;
p0      = 3.5*10^5; 
Sp_avg = 2 * Stroke * omega / 2 / pi;
pm = p0 * ((VDisp+Vc)/V)^gamma
w1 = C1*Sp_avg;
w2 = C2*(VDisp*T0)/(p0*V0)*(p-pm);

C=3.26;
macht=0.8;

alpha = C * Bore^(macht-1) * (p/1000)^macht * (w1+w2)^macht * T^(0.75-1.62*macht)

dQhl        = alfa*(A_c*(Twall-T)+A_p*(Twall-T));

%% DAE formulation
Rg = StateCyl.Rg;
yp = [dQhl-p*dVdt+hpaI*dmdtI+hpaE*dmdtE;
    p*V-Rg*T*m;
    dmidt];
end

