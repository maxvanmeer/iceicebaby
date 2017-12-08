clc
clear
close all

sheet1 = xlsread('WbData.xlsx', 'Settings');

sheet2 = xlsread('WbData.xlsx', 'WiebePar');

% each of the 10 cases has 6 rows

n=-1;
for i = 1:10 % loop over IDs
    allCases(i).CaseID = sheet1(n+2,1);
    
    allCases(i).RPM_sel = sheet1(n+2,3);
    allCases(i).RPM_act = sheet1(n+3,3);
    allCases(i).TORQUE = sheet1(n+4,3);
    
    allCases(i).T_plenum = sheet1(n+2,7);
    allCases(i).p_plenum = sheet1(n+3,7);
    allCases(i).lambda_AF = sheet1(n+4,7);
    allCases(i).lambda_exh = sheet1(n+5,7);
    allCases(i).EGRf = sheet1(n+6,7);
    n=n+6;
    
end

n = -1;
for i = 1:10
    allCases(i).f = [sheet2(n+2,7), sheet2(n+2,8), sheet2(n+2,9)];
    allCases(i).a = [sheet2(n+3,7), sheet2(n+3,8), sheet2(n+3,9)];
    allCases(i).n = [sheet2(n+3,7), sheet2(n+4,8), sheet2(n+4,9)];
    allCases(i).dCA = [sheet2(n+4,7), sheet2(n+5,8), sheet2(n+5,9)];
    allCases(i).CAign = [sheet2(n+5,7), sheet2(n+6,8), sheet2(n+6,9)];
    
    n = n+6;


end

% save(char('cases.mat'),a);
save('cases.mat','allCases');