clear all;
close all;
clc
p_all = 10:30;
for i = 1:length(p_all)
    p =    p_all(i);
    Temp = 850:1050;
    EGRf = 0.2;
    AP=0.405;
    nP=0.623;
    TaP=2977;
    mP=0.173;
    IDtP = AP*p^(-nP)*exp(TaP./Temp)*EGRf^mP;
    hold on
    plot(Temp,IDtP)
end
IDtlowcomp = AP*42.8^(-nP)*exp(TaP./781)*EGRf^mP
IDthighcomp = AP*58.8^(-nP)*exp(TaP./975)*EGRf^mP
%compare the plot with slide 37, the lower the line, the higher the pressure.
%As can be seen TSOI needs to be between 850 and 1100. However, if doing so,
%p_plenum (for case 124 2.238) is too low since IDt goes in a range of 2.5-5.5 ms.
%The range of IDt in the plot on slide 38 is 0.2-1.4 ms approximately.
%This means the pressure needs to be at least 10 times bigger.
%According to heyewood, page 544, the pressure should be between 10
%and 30 bar and the temperature between 770 K and 980 K. The constants found by
%Bart Somers equal the numbers of Stringer et al. the most. Looking at the
%table on page 545, other numbers can be found as well as a reference to,
%but I could not access it. The doc is called: 'Relative Roles of Premixed
%and Diffusion Burning in Diesel Combustion'. Since we have a common rail
%injector, which is an direct injection (DI) engine, two options remain.
%Our found relation is not depended on the speed. The temperature values of the
%'low compression ratio' are lower than the values named earlier. This condition results in
% 1.34 ms at our model, which does not match the value of Stringer (2.60 ms).
%However, if the conditions p = 58.8 and Temp = 975 K are used, the value
%of Stringer is 0.508 and ours is 0.5132. Thus these values are chosen and
% it is indeed correct that we have a high compression ratio.