clc
clear
close all

sheet1 = xlsread('WbData.xlsx', 'Settings')

sheet2 = xlsread('WbData.xlsx', 'WiebePar');

% each of the 10 cases has 6 rows

n=-1;
for i = 1:10 % loop over IDs
    a(i).CaseID = sheet1(n+2,1);
    
    a(i).RPM_sel = sheet1(n+2,3);
    a(i).RPM_act = sheet1(n+3,3);
    a(i).TORQUE = sheet1(n+4,3);
    
    a(i).T_plenum = sheet1(n+2,7);
    a(i).p_plenum = sheet1(n+3,7);
    a(i).lambda_AF = sheet1(n+4,7);
    a(i).lambda_exh = sheet1(n+5,7);
    a(i).EGRf = sheet1(n+6,7);
    n=n+6;
    
end

n = -1;
for i = 1:10
    a(i).f = [sheet2(n+2,7), sheet2(n+2,8), sheet2(n+2,9)];
    a(i).f = [sheet2(n+3,7), sheet2(n+3,8), sheet2(n+3,9)];
    a(i).f = [sheet2(n+3,7), sheet2(n+4,8), sheet2(n+4,9)];
    a(i).f = [sheet2(n+4,7), sheet2(n+5,8), sheet2(n+5,9)];
    a(i).f = [sheet2(n+5,7), sheet2(n+6,8), sheet2(n+6,9)];
    
    n = n+6;


end
