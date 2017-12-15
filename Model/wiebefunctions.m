function [HRR] = wiebefunctions(Ca)
CA = Ca;
%CA = -360:360;

%Per phase
f = [0.010 0.9  0.090];
a = [1.862 3.225 2.253];
n = [3 2 3];
dCA = [1.137 33.465 70.028];
CAI = [-1.954 -1.542 1.9];
Premix = SingleWiebe(CA,a(1),n(1),dCA(1),CAI(1));
Mixing = SingleWiebe(CA,a(2),n(2),dCA(2),CAI(2));
Late = SingleWiebe(CA,a(3),n(3),dCA(3),CAI(3));

%Total
HRR = f(1)*Premix + f(2)*Mixing + f(3)*Late;
end