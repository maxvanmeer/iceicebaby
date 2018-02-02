clc
clear
close all

N_sim = 15; % Square this for total amount

RPM = linspace(550,2150,N_sim);
T = linspace(100,2599,N_sim);

line1 = 1100+((2600-1400)/400).*(RPM-700);
line2=1e10;line3=1e10;
wrongIndices = find(T>line1);

omega = 2*pi*RPM/60;

%% Run parametrized cases
n=1;
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