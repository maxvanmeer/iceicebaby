clc
clear
close all

T = [300 300 1500 2500 2000];
w = [850 1000 1000 1000 1500]./60*2*pi;

for i = 4:length(T)
    iCase = 900+i;
    try
        MainDAE(T(i),w(i),iCase);
    catch
        disp(["Error occured on case ",num2str(i)," with T=",num2str(T(i))," w=",num2str(w(i))]);
    end
end