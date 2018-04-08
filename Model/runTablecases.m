clc
clear
close all

T = [300 300 1500 2500 2000];
w = [850 1000 1000 1000 1500]./60*2*pi;

for i = 1:length(T)
    iCase = 900+i;
    try
        MainDAE(T(i),w(i),iCase);
    catch e
        disp(['Error occured on case ',num2str(i),' with T=',num2str(T(i)),' w=',num2str(w(i))]);
        fprintf(1,'There was an error! The message was:\n%s',e.message);

    end
end