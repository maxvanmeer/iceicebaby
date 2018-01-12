function [HRR] = wiebefunctions(Ca)
global iCase
CA = Ca;
%CA = 0:360;
%Load parameters
iCase = 131;
CaseName = ['parameters' num2str(iCase,'%3.3i') '.mat'];
load('cases.mat');
%Per phase
f = allCases(iCase-121).f;
a = allCases(iCase-121).a;
n = allCases(iCase-121).n;
dCA = allCases(iCase-121).dCA;
CAign = allCases(iCase-121).CAign;
Premix = SingleWiebe(CA,a(1),n(1),dCA(1),CAign(1));
Mixing = SingleWiebe(CA,a(2),n(2),dCA(2),CAign(2));
Late = SingleWiebe(CA,a(3),n(3),dCA(3),CAign(3));

%Total
HRR = f(1)*Premix + f(2)*Mixing + f(3)*Late;
HR = trapz(HRR);
plot(CA,HRR);
hold on
end

