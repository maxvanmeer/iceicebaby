clc
clear



load('currentCase.mat');
currentCase.mode = 'couple';
iCase = 122;
currentCase.T = 2637.3;
currentCase.w = 1200/60*2*pi;
currentCase.EGRf = 0.2;
save('currentCase.mat','currentCase');

CA = -360:0.25:360;
PSOI = 51.7;
TSOI = 873.9;
for i = 1:length(CA)
   HRR_param(i) = wiebefunctions(CA(i),TSOI,PSOI);
end


load('cases.mat');
currentCase = allCases(iCase-121);
currentCase.mode = 'case';
save('currentCase.mat','currentCase');

for i = 1:length(CA)
   HRR_case(i) = wiebefunctions(CA(i),TSOI,PSOI);
end

figure
plot(CA,HRR_param)
hold on
plot(CA,HRR_case)
legend('Parametrized','Case')