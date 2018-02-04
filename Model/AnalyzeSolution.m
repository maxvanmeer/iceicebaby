clear all;clc; %close all;
%%
%% Add path to general functions and set Runiv
addpath('General');
global Runiv SpS
Runiv = 8.3144598;
%% Datadir just to show how you can organize output
DataDir = 'output'; 
%% Units, for convenience only
g=1e-3;
ms=1e-3; 
bara=1e5;
mm=1e-3;cm=1e-2;dm=0.1;
J = 1/(3.6e6);
liter = dm^3;
%% Select a loadcase
iCase = 122;        % Select a case to be analyzed
CaseName = ['Case' num2str(iCase,'%3.3i') '.mat'];
% CaseName = 'paramCase101.mat'; % For debugging
SaveName = fullfile(DataDir,CaseName);
load(SaveName);
fprintf('Read solution of Case %3i from %s\n',iCase,SaveName);
whos
%% Do your stuff
iSpSel = [3 4 6 7]; % Skip N2
t=time;p = y(:,1);T=y(:,2);mi=y(:,iSpSel);

RPM = Settings.N;

REVS = RPM/60;trev = 1/REVS;nREVS = (t(end)-t(1))/trev; 
it = find(t > (nREVS-2)*trev & t <= nREVS*trev);
% it = find(t > (0)*trev & t <= 2*trev);

% indices for all cycles in the simulation
i_per_cycle = length(t)/(nREVS/2);
 for i = 1:nREVS/2
     iti = ceil((i-1)*i_per_cycle+1e-1:1:i*i_per_cycle);
     it_all(:,i) = iti;
 end
 
  for i=1:size(it_all,2)
    ittt = it_all(:,i);
    [valueA, i_compA] = min(abs(V(ittt)));
    ppA = p(ittt);
    allMax(i) = max(ppA);
 end
 
% it = it_all(:,end);
 [value, i_comp] = min(abs(V(it)));


tp = t(it);
pp = p(it);
Tp = T(it);
mip = mi(it,:);

tp_all = t(it_all);
pp_all = p(it_all);

%% Pressure, Temperature and mass of component diagram
figure(1)
subplot(2,2,1)
plot(tp/ms,pp/bara);
xlabel('t [ms]');ylabel('p [bara]');
subplot(2,2,2)
plot(tp/ms,Tp);
xlabel('t [ms]');ylabel('T [K]');
subplot(2,2,[3 4])
plot(tp/ms,mip/g);
xlabel('t [ms]');ylabel('m [g]');
legend(yNames{iSpSel});

%% pV diagram
figure(2)
Vp = V(it);
Vp_all = V(it_all);
VDisp = max(V) - min(V);

i=1;
while p(i)<= p(i+floor(length(it)/2)) % Algorithm to find first intersection point
    p_intersect1 = p(i);
    V_intersect1 = V(i);
    i_comp = i; % index at which the compression starts; crossover point
    
    i=i+1;
end
% Assuming a sine in the volume, the next intersection of the pv diagram is at the same distance
i_exp = floor(length(it)/2) + i_comp; % index at which the expansion starts; crossover point 2
V_intersect2 = Vp(i_exp);
p_intersect2 = pp(i_exp);

Vcomp_exp = Vp(i_comp:i_exp);
pcomp_exp = pp(i_comp:i_exp);
V_pumploop = Vp([1:i_comp, i_exp:length(Vp)]);
p_pumploop = pp([1:i_comp, i_exp:length(Vp)]);

subplot(1,2,1)
hold on
plot(V_intersect1/liter,p_intersect1/bara,'g*') 
plot(V_intersect2/liter,p_intersect2/bara,'g*') 

pl=plot(V/liter,p/bara,'-',Vp/liter,pp/bara,'r-');
set(pl(end),'LineWidth',2);
xlabel('V [l]');ylabel('p [bara]');
subplot(1,2,2)
pl=loglog(V/liter,p/bara,'-',Vp/liter,pp/bara,'r-');
set(pl(end),'LineWidth',2);
xlabel('log V [l]');ylabel('log p [bara]');
%% Computations
W   = trapz(Vp,pp);     % Work, integral pdV
Wcomp_exp = trapz(Vcomp_exp,pcomp_exp);
W_pumploop = trapz(V_pumploop, p_pumploop);


for i = 1:size(it_all,2)
    W_all(1,i) = trapz(Vp_all(:,i),pp_all(:,i));
end

dummy = find(t > (nREVS-1.25)*trev); % I guess this is after IVC and before combustion. There are better ways.
index = dummy(1);
mfuel = mi(index,1);% fuel mass after intake valve close (just before combustion starts for instance)


figure(1)
subplot(2,2,[3 4])
line(t(index)*[1 1]/ms,mfuel*[1 1]/g,'Marker','o','MarkerSize',8,'MarkerFaceColor','y');
tx=text(t(index)*[1 1]/ms,1.1*mfuel*[1 1]/g,'Selected fuel mass','Rotation',45);
QLHV = Comb.QLHV;
Qin = mfuel*QLHV;


%% Torque for 6 cyclinders
T_all = W_all/(2*pi*(nREVS/Settings.Ncyc));     % Torque of every cycle
T_mean = sum(T_all)/(Settings.Ncyc);            % Mean torque of multiple cycles
T_V6 = 6*T_mean;                                % Mean torque of 6 cylinder engine
T_input = 100;
% mfuel2 = (2*2*pi*T_input/(0.46*QLHV)+0.00011);
% Qin2=mfuel2*QLHV;
% eff2=W/Qin2;

bsfc = mfuel/W*1000/J;
IMEP_net = W/VDisp;                                 % By work of second revolution of the cycle
IMEP_gross = Wcomp_exp/VDisp;                       % By work of complete cycle
PMEP = -W_pumploop/VDisp;                           % By work of first revolution (= IMEP_gross - IMEP_net)
BMEP = (2*pi*(nREVS/Settings.Ncyc)*T_mean)/VDisp;   % By work of measured brake torque
FMEP = IMEP_gross - BMEP;                           % By work caused by friction (= IMEP_gross - BMEP)

eff = W/Qin;
eff_all = sum(W_all)/Settings.Ncyc/Qin;


%% Plotting heat release rate
HRR_plot = HRR(350:420);
CA = -360:360;
CA_plot = 1:length(HRR_plot);

figure(3)
plot(CA,HRR)
xlabel('Crankangle [Degrees]')
ylabel('Heat release rate [J/degree]')
title('Heat release rate')
hold on
plot(CA10,HRR(361+CA10),'r*');
plot(CA50,HRR(361+CA50),'r*');
plot(CA90,HRR(361+CA90),'r*');
