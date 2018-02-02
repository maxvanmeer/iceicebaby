clc
clear
close all

% Open .mat files and make plots

load('paramInput.mat');




line1 = 1100+((2600-1400)/400).*(RPM_input-700);
line2=2600;
line3=1350+((1350-2850)/1000).*(RPM_input-2200);
line4=1050+((1050-2575)/600).*(RPM_input-2200);
wrongIndices = find(T_input>line1 | T_input>line2 | T_input>line3 | T_input>line4);

RPM_input(wrongIndices)=[];
T_input(wrongIndices)=[];

n=length(T_input);


for i = 1:n
    myCase(i) = load(['output/paramCase',num2str(i,'%3.3i'),'.mat']);
    t = myCase(i).time;
    REVS = myCase(i).Settings.N/60;
    trev = 1/REVS;
    nREVS = (t(end)-t(1))/trev;
    it = find(t > (nREVS-2)*trev & t <= nREVS*trev);
    
    p = myCase(i).y(:,1);
    
    Vp = myCase(i).V(it);
    pp = p(it);
    
    W   = trapz(Vp,pp); % Work, integral pdV
    QLHV = myCase(i).Comb.QLHV;
    
    iSpSel = [3 4 6 7]; % Skip N2
    
    dummy = find(t > (nREVS-1.25)*trev); % I guess this is after IVC and before combustion. There are better ways.
index = dummy(1);

    mi=myCase(i).y(:,iSpSel);
    mfuel = mi(index,1);% fuel mass after intake valve close (just before combustion starts for instance)
    
    Qin = mfuel*QLHV;
    
    
    efficiency(i) = W/Qin;
    
end

dRPM = (max(RPM_input)-min(RPM_input))/1000;
dT = (max(T_input)-min(T_input))/1000;
RPM_range = 500:dRPM:2250;
T_range = 50:dT:2250;
F = scatteredInterpolant(RPM_input',T_input',efficiency','natural','none');
[RPMq,Tq] = meshgrid(RPM_range,T_range);
q = F(RPMq,Tq);

figure
mesh(RPMq,Tq,q);
xlabel('RPM [-]');
ylabel('Torque [Nm]');
title('Efficiency map');
cc=colorbar;
grid
set(cc, 'Fontsize', 20)
set(gca,'FontSize',20)