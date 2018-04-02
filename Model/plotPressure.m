clc

clear

%% set parameters
iCase = 127;
Torque = 1325.1;
w = 1200/60*2*pi;
%% run case
MainDAE(iCase);


%% Units, for convenience only
g=1e-3;
ms=1e-3; 
bara=1e5;
mm=1e-3;cm=1e-2;dm=0.1;
J = 1/(3.6e6);
liter = dm^3;
%% Do your stuff (case)

CaseName = ['Case',num2str(iCase),'.mat']; % For debugging

DataDir = 'output'; 
SaveName = fullfile(DataDir,CaseName);
load(SaveName);

iSpSel = [3 4 6 7]; % Skip N2
t=time;p = y(:,1);T=y(:,2);mi=y(:,iSpSel);

RPM = Settings.N;

REVS = RPM/60;trev = 1/REVS;nREVS = (t(end)-t(1))/trev; 
it = find(t > (nREVS-2)*trev & t <= nREVS*trev); % Last cycle
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


%% run Parametrized
MainDAE(Torque,w,iCase);
clearvars -except iCase DataDir tp pp
%% Units, for convenience only
g=1e-3;
ms=1e-3; 
bara=1e5;
mm=1e-3;cm=1e-2;dm=0.1;
J = 1/(3.6e6);
liter = dm^3;
%% Do your stuff (case)

CaseName = ['paramCase',num2str(iCase),'.mat']; % For debugging

SaveName = fullfile(DataDir,CaseName);
load(SaveName);

iSpSel = [3 4 6 7]; % Skip N2
t=time;p = y(:,1);T=y(:,2);mi=y(:,iSpSel);

RPM = Settings.N;

REVS = RPM/60;trev = 1/REVS;nREVS = (t(end)-t(1))/trev; 
it = find(t > (nREVS-2)*trev & t <= nREVS*trev); % Last cycle
% it = find(t > (0)*trev & t <= 2*trev);

% indices for all cycles in the simulation
i_per_cycle = length(t)/(nREVS/2);
 for i = 1:nREVS/2
     iti = ceil((i-1)*i_per_cycle+1e-1:1:i*i_per_cycle);
     it_all(:,i) = iti;
 end
 

 
% it = it_all(:,end);
 [value, i_comp] = min(abs(V(it)));


tp2 = t(it);
pp2 = p(it);
Tp = T(it);
mip = mi(it,:);



%% Plot both
figure
hold on
% plot((tp-tp(1))/ms,pp/bara);
% plot((tp2-tp(1))/ms,pp2/bara);
plot(linspace(-360,360,length(pp)),pp/bara);
plot(linspace(-360,360,length(pp)),pp2/bara);
title(['Pressure trace of case ',num2str(iCase)]);
legend('Case','Parametrized');
xlim([-360 360]);
ylabel('p [bar]');
xlabel('Crank angle [deg]');
grid


