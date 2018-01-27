clc
clear
close all

RPM = linspace(500,2250,5);
T = linspace(50,2250,4);

omega = 2*pi*RPM/60;

% Run parametrized cases
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