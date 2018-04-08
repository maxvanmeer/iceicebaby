function [HRR] = wiebefunctions(Ca,TSOI,PSOI)
CA = Ca;
%CA = -360:360;

load('currentCase.mat');
mode = currentCase.mode;

%Per phase
global CAignP EOId

if strcmp(mode,'case')
    % Do wiebe with case data
    f = currentCase.f;
    a = currentCase.a;
    nP = currentCase.n;
    dCA = currentCase.dCA;
    CAI = currentCase.CAign;
    
    Premix = SingleWiebe(CA,a(1),nP(1),dCA(1),CAI(1));
    Mixing = SingleWiebe(CA,a(2),nP(2),dCA(2),CAI(2));
    Late = SingleWiebe(CA,a(3),nP(3),dCA(3),CAI(3));
    HRR = f(1)*Premix + f(2)*Mixing + f(3)*Late;
    
elseif strcmp(mode,'couple')
    % Do wiebe with parametrization
    T = currentCase.T;
    w = currentCase.w;
    EGRf= currentCase.EGRf;
    EOIt = deg2rad(-2+22*T*6/2600)/w*1000;    %[ms]!!
    INJ_durt = 0.5+3*T*6/2600;                %[ms]
    SOIt = EOIt-INJ_durt;                   %[ms]
    INJ_durd = rad2deg((INJ_durt/1000)*w);  %[CAD]
    SOId = rad2deg((SOIt/1000)*w);          %[CAD]
    EOId = rad2deg((EOIt/1000)*w);          %[CAD]
    % Parameters needed for dP = pressure common rail? Can also be used to
    % calculate QLHV if dP is assumed to be a certain value.
    QLHV = 4.26e7;

    mfuel = (2*2*pi*T*6/(0.46*QLHV)+0.00011)/6;               % /6, values are for 6 cylinders. from hints
    
    Cd = 0.8;                                           %[-]
    D_hole = 180e-6;                                    %[m]
    A_holes = 7*(pi/4)*D_hole^2;                        %[m2]
    rho_D = 0.810*10^3;                                 %[kg/m3]
    m_rate = (mfuel)/(INJ_durt/1000);                   %[kg/s]
    dP = ((m_rate)/(A_holes*Cd*sqrt(rho_D)))^2*10^-5;   %[bar]
    ID_alt = (0.35 + 0.5*(sqrt(dP)-sqrt(2500))/25)^2;   %[ms]
    %ID_alt is calculated by the formula used in the appendix of the
    %studyguide. This method could be out of date but it does give a good
    %indication.
    
    %Argumentation for using these values can be found in 'Proof_TSOI_PSOI'
    %     Temp = 975;                         %[K]
    %     p = 58.8;                           %[bar]
    Temp = TSOI;
%     Temp = 880;
    p = PSOI;
%     p = 46;

%     Temp = 975;
%     p = 58.8;
    
    %Premix
    AP=0.405;
    nP=0.623;
    TaP=2977;
    mP=0.173;
    IDtP = AP*p^(-nP)*exp(TaP/Temp)*EGRf^mP; %[ms] Ignition Delay
%     IDtP = ID_alt;
    rt = IDtP / INJ_durt;   %[-]
    CAignP = rad2deg((SOIt+IDtP)/1000*w); % [deg=CAD]
    
    %Mix (same values as premix)
    CAignM = CAignP;
    
    %Late
    AL=0.237;
    nL=0.379;
    TaL=3289;
    mL=0.145;
    IDtL = AL*p^(-nL)*exp(TaL/Temp)*EGRf^mL; %[ms]
    CAignL = rad2deg(w*(IDtL+SOIt)/1000); % [CAD]
    
    %Calculate fP and fM
    if rt<0.8
        fP=-0.0107+0.1684*rt;
        fM=0.8793-0.3428*rt;
    else
        fP=-0.566+0.7627*rt;
        fM=1.0125-0.4228*rt;
    end
    
    %Calculate and clip fL
    fL = 1-fP-fM;
    if fL<0
        fP = fP-abs(fL)/2;
        fM = fM-abs(fL)/2;
        fL=0;
    end
    
    %Calculate dCA
    A_dCA=0.00405;
    n_dCA=0.631;
    Ta_dCA=6353.9;
    m_dCA=0.216;
    
    dtP = A_dCA*p^(-n_dCA)*exp(Ta_dCA/Temp)*EGRf^m_dCA;
    dCAP = rad2deg(dtP/1000*w); % [CAD]
    dCAM = -0.024*INJ_durd^2+1.1776*INJ_durd+10.7056; %[CAD]
    dCAL = -0.0473*INJ_durd^2+2.105*INJ_durd+21.657; %[CAD]
    
    %Calculate a, n
    aP = -0.534 *IDtP+1.878;
    nP_wiebe=3;
    
    aM=0.0031*INJ_durd^2-0.226*INJ_durd+4.901;
    nM_wiebe=2;
    
    lna = -0.135*INJ_durd+2.25;
    aL = exp(lna);
    nL_wiebe=3;
    
    Premix = SingleWiebe(CA,aP,nP_wiebe,dCAP,CAignP);
    Mixing = SingleWiebe(CA,aM,nM_wiebe,dCAM,CAignM);
    Late = SingleWiebe(CA,aL,nL_wiebe,dCAL,CAignL);
    HRR = fP*Premix + fM*Mixing + fL*Late;
end


end