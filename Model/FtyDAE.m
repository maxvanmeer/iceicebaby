function [yp] = FtyDAE(t,y)
global Int Exh QLHV SpS Runiv omega Qheatloss dt
global  mfuIVCClose si EtaComb Bore Stroke rc CA50 BDUR SOC EOC CAignP EOId

Twall   = 273+80;   %80 degrees Celcius is on the lower side. Higher is better for the engine efficieny.
Tpiston = 273+110;  %110 degrees Celsius is a guess, based on the normal temperature of engine oil (which cools the pistons). NOT SURE.
% alfa    = 500;

VDisp   = pi*(Bore/2)^2*Stroke;
Vc      = VDisp/(rc-1);

load('currentCase.mat');


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
Ca          = time2ang(t,omega)/(2*pi)*360; %Crankangle vector for the entire simulation time
reducedCa   = mod(Ca+360,720)-360;          %Crankangle for one cycle
CADS        = omega/(2*pi)*360;             %Rotational speed in degrees per second
%% Taking parameters from NASA database
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
% dmfuComb    = EtaComb*mfuel*wiebefunctions(reducedCa)*CADS;


global Torque w TSOI PSOI haveToSetSOI
EOIt = deg2rad(-2+22*Torque/2600)/w*1000;    %[ms]!!
INJ_durt = 0.5+3*Torque/2600;                %[ms]
SOIt = EOIt-INJ_durt;                   %[ms]
SOId = rad2deg((SOIt/1000)*w);          %[CAD]
EOId = rad2deg((EOIt/1000)*w);          %[CAD]


if strcmp(currentCase.mode,'couple')
    if reducedCa <= -350
        haveToSetSOI = true;
    end
    if reducedCa >= SOId && haveToSetSOI == true% You're after start of injection
        TSOI = T;
        PSOI = p*1e-5;
        haveToSetSOI = false;
    end
end



dmfuComb    = EtaComb*mfuIVCClose*wiebefunctions(reducedCa,TSOI,PSOI)*CADS;


dmidt_c     = si'*dmfuComb;
dmidt       = dmidt - dmidt_c;
dQcomb      = QLHV*dmfuComb;
dQcomb_real = ei*dmidt_c;

%% Woschni heatloss
%normal model
% CA50=10;                        % CA50 (50% HR), replace with
HR = currentCase.HR;
ReducedCA = -360:360;

% Useful and important Crank Angles determined from heat release rate:
CA01 = ReducedCA(find(HR>0.01*HR(length(HR)),1)+1); 
CA10 = ReducedCA(find(HR>0.1*HR(length(HR)),1)+1);
CA50 = ReducedCA(find(HR>0.5*HR(length(HR)),1)+1);
CA90 = ReducedCA(find(HR>0.9*HR(length(HR)),1)+1);
% BDUR=20;                        % Burn Duration, replace with data case
BDUR = CA90-CA10;
CAign = CA50-0.5*BDUR;
CAend = CAign+BDUR;

if haveToSetSOI
    SOC = CAign;
    EOC = CAend;
else
    SOC = CAignP;
    EOC = EOId;
end


Tr      = T_plenum;             % Plenum temperature
pr      = p_plenum;             % Plenum pressure

% V0      = CylVolumeFie(t(1));
% T0      = 273;
% p0      = 3.5*10^5;
global T_plenum p_plenum LCon Vmax

pm = ((VDisp+Vc)/(V))^gamma * p0;

alfa = alfaWoschni(reducedCa,T,p,pm,Tr,pr,Vr);

alfaplot = [alfaplot alfa];

dQhl= alfa*(A_c*(Twall-T)+A_p*(Tpiston-T));

Qheatloss = Qheatloss + dQhl * dt;

%% DAE formulation
Rg = StateCyl.Rg;
yp = [dQhl-p*dVdt+hpaI*dmdtI+hpaE*dmdtE;
    p*V-Rg*T*m;
    dmidt];
end
