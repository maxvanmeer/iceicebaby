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
reducedCa   = mod(Ca+360,720)-360
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
%normal model
CA50=10;                        % CA50 (50% HR), replace with 
BDUR=20;                        % Burn Duration, replace with data case
CAign = CA50-0.5*BDUR;
CAend = CAign+BDUR;

%for the intake/exhaust Ca, look at 'valvedata.mat'
if reducedCa >= CAign && reducedCa <= CAend         %combustion phase
    C1 = 2.28;
    C2 = 3.24*10-3;
elseif (reducedCa >= -360 && reducedCa <= -108) || (reducedCa >= 91 && reducedCa <= 360)      %Intake/Exhaust phase
    C1 = 6.18;
    C2 = 0;
else
    C1 = 2.28;
    C2 = 0;
end

V0      = CylVolumeFie(t(1));
T0      = 273;
p0      = 3.5*10^5;                                             %might need to be adjusted because of looad change
Sp_avg = 2 * Stroke * omega / 2 / pi;
ratio = 1 + VDisp/Vc;
pm = ((VDisp+Vc)/(V))^gamma * p0;      
w = C1 * Sp_avg + C2 * (VDisp*T0)/(p0*V0) * (p-pm)/1000;       
if (p-pm)<0
    w = 0;
end
alpha=3.26;
macht=0.8;
hc_tot = alpha * Bore^(macht-1) * (p/1000)^macht * (w)^macht * T^(0.75-1.62*macht)


% %Heat loss HCCI engines, due to new paper given by Somers
% macht_HCCI = 0.8
% alpha_HCCI = 3.26 %? Or 0.34?
% 
% w_HCCI = C1 * Sp_avg + C2/6 * (VDisp*T0)/(p0*V0) * (p-pm)/1000    %adjusted woschni because new HCCI paper canvas, again unsure about pm
% if (p-pm2)<0
%     w_HCCI = 0;
% end
% hc_tot_HCCI = alpha_HCCI * ? * (p/1000)^macht_HCCI * (w2)^macht_HCCI * T^(0.73)          %adjusted alpha because new HCCI paper canvas
% 
% 
% dQhl        = hc_tot*(A_c*(Twall-T)+A_p*(Twall-T));

%% DAE formulation
Rg = StateCyl.Rg;
yp = [dQhl-p*dVdt+hpaI*dmdtI+hpaE*dmdtE;
    p*V-Rg*T*m;
    dmidt];
end

