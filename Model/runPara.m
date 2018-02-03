clc
clear
close all

N_sim = 4; % Square this for total amount

RPM = linspace(550,2150,N_sim);
T = linspace(100,2599,N_sim);

line1 = 1100+((2600-1400)/400).*(RPM-700);
line2=2600;
line3=1350+((1350-2850)/1000).*(RPM-2200);
line4=1050+((1050-2575)/600).*(RPM-2200);
wrongIndices = find(T>line1 | T>line2 | T>line3 | T>line4);

% RPM = RPM(setdiff(1:end,wrongIndices)); % Strip wrong indices
% T = T(setdiff(1:end,wrongIndices));

omega = 2*pi*RPM/60;



%% Run parametrized cases
n=1;
RPM_input=[];
T_input=[];
for i=1:length(omega)
    for j=1:length(T)
        MainDAE(T(j),omega(i),n);
        RPM_input(n) = RPM(i);
        T_input(n) = T(j);
        
        n=n+1;
    end
end
save('paramInput.mat','RPM_input','T_input');

postProcess