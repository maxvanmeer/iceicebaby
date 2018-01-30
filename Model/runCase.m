clc
clear
close all

global iCase 
iCase = input('Select your loadcase (122 t/m 131). Enter 0 to run all cases: ');


if (iCase==0) 
  for i = 122:131;
  MainDAE(i); 
   end
else
   MainDAE(iCase); 
end