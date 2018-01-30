function [HRR] = wiebefunctions(Ca,myCase)

    load('cases.mat');
    % Do wiebe with case data
    f = allCases(myCase-121).f;
    a = allCases(myCase-121).a;
    n = allCases(myCase-121).n;
    dCA = allCases(myCase-121).dCA;
    CAign = allCases(myCase-121).CAign;
    
    Premix = SingleWiebe(Ca,a(1),n(1),dCA(1),CAign(1));
    Mixing = SingleWiebe(Ca,a(2),n(2),dCA(2),CAign(2));
    Late = SingleWiebe(Ca,a(3),n(3),dCA(3),CAign(3));
    HRR = f(1)*Premix + f(2)*Mixing + f(3)*Late;

end

