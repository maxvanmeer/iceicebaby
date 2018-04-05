function [M] = MassDAE( t,y )
% % y
% if isnan(y(1))
%     
% 
%     M=[0 0.504545503677356 -1202484.61492994 -63045.0982699346 -72230.5969746713 -8983101.99482530 -13530139.5474623;0 0 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1]
%     return;
% end
global SpS
% Mass matrix belonging to MainDAE
%   Determines M of M dydt= Fty;
%   Assumes : y = [p Temp m1 ... m5]; Fty = [EnConv; Eq of State;dm1dt;...;dm5dt];
p=y(1);Temp=y(2);mi=y(3:end); %m=sum(mi);
V=CylVolumeFie(t);
Nsp = length(SpS);
Cvi     = mi; % Will be overwritten
ei      = mi; % idem
for ii=1:Nsp
    Cvi(ii) = CvNasa(Temp,SpS(ii));
    ei(ii) = ENasa(Temp,SpS(ii));
end
mCv = Cvi'*mi;           % m * Cv
%%
M = [0 mCv ei';...
    0   0  0 0 0 0 0;...
    0   0  1 0 0 0 0;...
    0   0  0 1 0 0 0;...
    0   0  0 0 1 0 0;...
    0   0  0 0 0 1 0;...
    0   0  0 0 0 0 1;...
    ];
end

