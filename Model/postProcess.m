clc
clear
close all

% Open .mat files and make plots

load('paramInput.mat');

allParamCases = 1:length(T_input);


line1 = 1100+((2600-1400)/400).*(RPM_input-700);
line2=2600;
line3=1350+((1350-2850)/1000).*(RPM_input-2200);
line4=1050+((1050-2575)/600).*(RPM_input-2200);
wrongIndices = find(T_input>line1 | T_input>line2 | T_input>line3 | T_input>line4);

RPM_input(wrongIndices)=[];
T_input(wrongIndices)=[];

goodParams = setdiff(allParamCases,wrongIndices);


J = 1/(3.6e6);


for i =1:length(goodParams)
    myParam = goodParams(i);
    myCase(i) = load(['output/paramCase',num2str(myParam,'%3.3i'),'.mat']);
    t = myCase(i).time;
    V = myCase(i).V;
    REVS = myCase(i).Settings.N/60;
    trev = 1/REVS;
    nREVS = (t(end)-t(1))/trev;
    it = find(t > (nREVS-2)*trev & t <= nREVS*trev);
    
    % indices for all cycles in the simulation
    i_per_cycle = length(t)/(nREVS/2);
    for j = 1:nREVS/2
        iti = ceil((j-1)*i_per_cycle+1e-1:1:j*i_per_cycle);
        it_all(:,j) = iti;
    end
    
    p = myCase(i).y(:,1);
    
    Vp = myCase(i).V(it);
    Vp_all = myCase(i).V(it_all);
    pp = p(it);
    
     tp_all = t(it_all);
 pp_all = p(it_all);
    
    j=1;
    while p(j)<= p(j+floor(length(it)/2)) % Algorithm to find first intersection point
        p_intersect1 = p(j);
        V_intersect1 = V(j);
        i_comp = j; % index at which the compression starts; crossover point
        
        j=j+1;
    end
    % Assuming a sine in the volume, the next intersection of the pv diagram is at the same distance
    i_exp = floor(length(it)/2) + i_comp; % index at which the expansion starts; crossover point 2
    V_intersect2 = Vp(i_exp);
    p_intersect2 = pp(i_exp);
    
    Vcomp_exp = Vp(i_comp:i_exp);
    pcomp_exp = pp(i_comp:i_exp);
    V_pumploop = Vp([1:i_comp, i_exp:length(Vp)]);
    p_pumploop = pp([1:i_comp, i_exp:length(Vp)]);
    
    
    
    W   = trapz(Vp,pp); % Work, integral pdV
    Wcomp_exp = trapz(Vcomp_exp,pcomp_exp);
    for j = 1:size(it_all,2)
        W_all(1,j) = trapz(Vp_all(:,j),pp_all(:,j));
    end
    
    QLHV = myCase(i).Comb.QLHV;
    
    iSpSel = [3 4 6 7]; % Skip N2
    
    dummy = find(t > (nREVS-1.25)*trev); % I guess this is after IVC and before combustion. There are better ways.
    index = dummy(1);
    
    mi=myCase(i).y(:,iSpSel);
    mfuel = mi(index,1);% fuel mass after intake valve close (just before combustion starts for instance)
    
    Qin = mfuel*QLHV;
    
    bsfc_real(i) = mfuel/W*1000/J;
    w = RPM_input(i)/60*2*pi;
    bsfc_fakenews(i) = mfuel/(T_input(i)*w*2*trev)*1000/J;
    
    efficiency(i) = W/Qin;
    
    T_all = W_all/(2*pi*(nREVS/myCase(i).Settings.Ncyc));     % Torque of every cycle
    T_mean = sum(T_all)/(myCase(i).Settings.Ncyc);            % Mean torque of multiple cycles
    
    myW(i) = W;
    VDisp = max(V) - min(V);
    IMEP_net(i) = W/VDisp;                                 % By work of second revolution of the cycle
    IMEP_gross(i) = Wcomp_exp/VDisp;                       % By work of complete cycle
    BMEP(i) = (2*pi*(nREVS/myCase(i).Settings.Ncyc)*T_mean)/VDisp;   % By work of measured brake torque
    FMEP(i) = IMEP_gross(i) - BMEP(i);                           % By work caused by friction (= IMEP_gross - BMEP)

    
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
xlabel('RPM [1/min]');
ylabel('Torque [Nm]');
title('Efficiency map');
cc=colorbar;
grid
% set(cc, 'Fontsize', 20)
% set(gca,'FontSize',20)
hold on
plot3(RPM_input,T_input,efficiency*1.5,'o');

%bsfc

F_bsfc = scatteredInterpolant(RPM_input',T_input',bsfc_real','natural','none');
q_bsfc = F_bsfc(RPMq,Tq);

figure
mesh(RPMq,Tq,q_bsfc);
xlabel('RPM [1/min]');
ylabel('Torque [Nm]');
title('bsfc map');
cc=colorbar;
grid
% set(cc, 'Fontsize', 20)
% set(gca,'FontSize',20)

%imep
F_imep = scatteredInterpolant(RPM_input',T_input',IMEP_gross','natural','none');
q_imep = F_imep(RPMq,Tq);


figure
mesh(RPMq,Tq,q_imep);
xlabel('RPM [1/min]');
ylabel('Torque [Nm]');
title('imep');
cc=colorbar;
grid
% set(cc, 'Fontsize', 20)
% set(gca,'FontSize',20)

%bmep
F_BMEP = scatteredInterpolant(RPM_input',T_input',BMEP','natural','none');
q_BMEP = F_BMEP(RPMq,Tq);

figure
mesh(RPMq,Tq,q_BMEP);
xlabel('RPM [1/min]');
ylabel('Torque [Nm]');
title('bmep');
cc=colorbar;
grid
% set(cc, 'Fontsize', 20)
% set(gca,'FontSize',20)

%fmep
F_FMEP = scatteredInterpolant(RPM_input',T_input',FMEP','natural','none');
q_FMEP = F_FMEP(RPMq,Tq);

figure
mesh(RPMq,Tq,q_FMEP);
xlabel('RPM [1/min]');
ylabel('Torque [Nm]');
title('fmep');
cc=colorbar;
grid
% set(cc, 'Fontsize', 20)
% set(gca,'FontSize',20)
