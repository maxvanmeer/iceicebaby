function [HRR] = wiebefunctions(Ca)
CA = Ca;
%CA = -360:360;

load('currentCase.mat');
mode = currentCase.mode;

%Per phase
% f = [0.010 0.9  0.090];
% a = [1.862 3.225 2.253];
% n = [3 2 3];
% dCA = [1.137 33.465 70.028];
% CAI = [-1.954 -1.542 1.9];

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
    EOI = -2+18*T/2700;
    INJ_dur = 0.5+3*T/2700;
    SOI = EOI-INJ_dur;
    
    Temp = 1050; %no argumentation
    p = currentCase.p_plenum;
    
    %Premix
    AP=0.405;
    nP=0.623;
    TaP=2977;
    mP=0.173;
    IDtP = AP*p^(-nP)*exp(TaP/Temp)*EGRf^mP;
    rtP = IDtP / INJ_dur;
    CAignP = w*IDtP/1000; % ASSUMPTION: no angle offset needed
    
    %Mix (apparently same values as premix)
    AM=0.405;
    nM=0.623;
    TaM=2977;
    mM=0.173;
    IDtM = AM*p^(-nM)*exp(TaM/Temp)*EGRf^mM;
    rtM = IDtM / INJ_dur;
    CAignM = w*IDtM/1000; % ASSUMPTION: no angle offset needed
    
    %Late
    AL=0.237;
    nL=0.379;
    TaL=3289;
    mL=0.145;
    IDtL = AL*p^(-nL)*exp(TaL/Temp)*EGRf^mL;
    rtL = IDtL / INJ_dur;
    CAignL = w*IDtL/1000; % ASSUMPTION: no angle offset needed

    
    
    %Calculate f
    if rtP<0.8
        fP=-0.0107+0.1684*rtP;
    else
        fP=-0.566+0.7627*rtP;
    end
    
    if rtM < 0.8
        fM=-0.0107+0.1684*rtM;
    else
        fM=1.0125-0.4228*rtP;
    end
    
    %Calculate and clip fL
    fL = 1-fP-fM;
    if fL<0
        fP = fP-abs(fL)/2;
        fM = fM-abs(fL)/2;
        fL=0;
    end
    
    %Calculate dCa
    A_dCa=0.00405;
    n_dCa=0.631;
    Ta_dCa=6353.9;
    m_dCa=0.216;
    
    dtP = A_dCa*p^(-n_dCa)*exp(Ta_dCa/Temp)*EGRf^m_dCa;
    dCaP = dtP/1000*w; % ASSUMPTION: no angle offset needed
    dCaM = -0.024*INJ_dur^2+1.1776*INJ_dur+10.7056;
    dCaL = -0.0473*INJ_dur^2+2.105*INJ_dur+21.657;
    
    %Calculate a, n
    aP = -0.534 *IDtP+1.878;
    nP_wiebe=3;
    
    aM=0.0031*INJ_dur^2-0.226*INJ_dur+4.901;
    nM_wiebe=2;
    
    aL = exp(-0.135*INJ_dur+2.25);
    nL_wiebe=3;
    
    Premix = SingleWiebe(CA,aP,nP_wiebe,dCaP,CAignP);
    Mixing = SingleWiebe(CA,aM,nM_wiebe,dCaM,CAignM);
    Late = SingleWiebe(CA,aL,nL_wiebe,dCaL,CAignL);
    HRR = fP*Premix + fM*Mixing + fL*Late;
    
    fs=[fP,fM,fL];
%     save('wiebe.mat','premix','mix','late','fs');
end


end