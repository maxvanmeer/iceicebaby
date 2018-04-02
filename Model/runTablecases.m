clc
clear
close all

T = [300 300 1500 2500 2000];
w = [850 1000 1000 1000 1500]./60*2*pi;

for i = 4:length(T)
    iCase = 900+i;
    MainDAE(T(i),w(i),iCase);
end