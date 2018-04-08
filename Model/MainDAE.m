function MainDAE(varargin)
global reducedCaPlot CaPlot
reducedCaPlot = [];
CaPlot = [];
%% MainDAE %%
% The MainDAE file is able to run both a single case, or a T-w couple.

%If you only want to run a case, type MainDAE(iCase) into the command window

%If you want to run a T-w couple, type MainDAE(T,w) into the command window
%1st arg: Input torque  T  in [Nm]
%2nd arg: Angular velocity  w  in [rad/s]

%Give only one argument if you want to use a case, or two arguments if you
%want to use a T-w couple

%Running the model manually runs the default case which is defined here:
mode = 'case';
load('cases.mat');
defaultCase = 122; % In case you want to run this file directly

global w Torque

if nargin == 0
    iCase = defaultCase;
    
elseif nargin == 1 % Run a specific case
    iCase = varargin{1};
    
elseif nargin ==3
    T = varargin{1}/6;
    w = varargin{2};
    Torque = T;

    paramNumber = varargin{3};
    mode = 'couple';
    
    
elseif nargin==2
    T=varargin{1}/6;
    w=varargin{2};
    Torque = T;

    paramNumber=1;
    mode = 'couple';
end

if (~exist('iCase','var'))
    iCase = defaultCase; %Catching exception
end

switch mode
    case 'case'
        disp(['Running case ',num2str(iCase),'...']);
        
    case 'couple'
        disp(['Running couple T=',num2str(T),' w=',num2str(w),'...']);
end


%% Add path to general functions and set Runiv
addpath('General');
global Runiv SpS QLHV p_plenum T_plenum Vr
global CA01 CA05 CA10 CA50 CA90 CA95 BDUR
Runiv = 8.3144598;
%% Datadir just to show how you can organize output
DataDir = 'output';
%% Engine Parameters
% global iCase
% iCase = input('Select your loadcase (122 t/m 131))');
CaseName = ['cases.mat'];
load(CaseName);

%% Units, for convenience only
g=1e-3;
ms=1e-3;
bara=1e5;
mm=1e-3;cm=1e-2;dm=0.1;
liter = dm^3;
%% Set a few global variables
global rc LCon Stroke Bore N omega Di De VDisp dt p_plenum T_plenum Vmax % Engine globals
LCon    = 261.6*mm;                 % connecting rod length
Stroke  = 158*mm;                   % stroke
Bore    = 130*mm;                   % bore
rc      = 17.45;                    % compression ratio
VDisp   = pi*(Bore/2)^2*Stroke;     % Displacement Volume

global TSOI PSOI haveToSetSOI
TSOI = -13.33*log(T) + 943.89
PSOI = 0.0368*T + 40.752


haveToSetSOI = true;

if strcmp(mode,'case')
    N       = allCases(iCase-121).RPM_act;          % N represents the RPM
elseif strcmp(mode,'couple')
    N = w/(2*pi)*60; %revs/s
    
end
Cyl.LCon = LCon;Cyl.Stroke=Stroke;Cyl.Bore=Bore;Cyl.rc=rc;

%% Simple combustion model settings (a gaussian distribution)
global mfuIVCClose si EtaComb Qheatloss alfaplot

CA50=10;                        % Gaussian value, no longer relevant
BDUR=20;                        % Gaussian value, no longer relevant
EtaComb = 0.99;                 % Combustion efficiency

%% Intake and exhaust pressures

if strcmp(mode,'case')
    p_plenum  = allCases(iCase-121).p_plenum*bara;       % plenum pressure
    T_plenum  = allCases(iCase-121).T_plenum;            % plenum temperature
    p_exhaust = p_plenum+0.1*bara;                       % exhaust back-pressure
elseif strcmp(mode,'couple')
    T_plenum = 320;
    %p_plenum will be calculated later on
    %p_exhaust too
end
T_exhaust  = 400;

%% Geometric and timing data of the valves
filenameValveData=fullfile('General','Valvedata.mat');
load(filenameValveData);
global Int Exh
Di = 48.5*mm;      % diameter inlet valve
De = 43.25*mm;     % diameter exhaust valve


%% Chemistry and Thermodynamic properties
filenameThermalDataBase=fullfile('General','NasaThermDatFull.mat');
load(filenameThermalDataBase);
indexes = myfind({Sp.Name},{'Diesel','O2','N2','CO2','H2O'});   % Means that the order is now set for all vectors if you want to use SpS lateron directly
SpS     = Sp(indexes);
Nsp     = length(SpS);
Names   = {SpS.Name};
Mi      = [SpS.Mass];
Xair    = [0 0.21 0.79 0 0];                                    % Air comp
Xfuel   = [1 0 0 0 0];                                          % Fuel comp
Yair    = Xair.*Mi/(Xair*Mi');
Yfuel   = [1 0 0 0 0];
nC      = SpS(1).Elcomp(3);                                     % SpS(1) is fuel by definition (line 50)
nH      = SpS(1).Elcomp(2);
nui     = [1   nC+nH/4 0 -nC -nH/2];                            % Reaction stoichiometry (molar)
si      = nui.*Mi/Mi(1);                                        % Reaction stoichiometry (mass)
AFstoi_molar  = nui(2)+nui(2)*Xair(3)/Xair(2);                  % So-called stoichiometric air fuel ratio (fuel property for given air composition), sometimes students use this as AFstoi.
AFstoi  = si(2)+si(2)*Yair(3)/Yair(2);                          % So-called stoichiometric air fuel ratio (fuel property for given air composition)
%% Set simulation time
Ncyc    = 20;
REVS    = N/60;
omega   = REVS*2*pi;
tcyc    = (2/REVS);
t       = [0:0.1:360]./360*tcyc*Ncyc;
Vmax    = max(CylVolumeFie(t));
dt      = t(100)-t(99);
CADS    = omega/(2*pi)*360;
%% Compute initial conditions and intake/exhaust composition
V0      = CylVolumeFie(t(1));
T0      = T_plenum;
% T0 = 273;

if strcmp(mode,'case')
    lambda  = allCases(iCase-121).lambda_AF;
    AF      = AFstoi*lambda;
    EGRf = allCases(iCase-121).EGRf/100;
elseif strcmp(mode,'couple')
    QLHV=4.26e7;
    mfuel = (2*2*pi*T*6/(0.46*QLHV)+0.00011)/6;               % /6, values are for 6 cylinders. from hints
    mair = (10.8e-3 + 28.4e-3 *T*6/2700)/6;                   % From hints, causes low p_plenum (26%)
%     mfuel = (2*pi*T*6/(2*0.46*QLHV) + 0.000027);                % NOT /6, mfuel and mair for 1 cylinder
%     mair = (2.7e-3 + 6.9e-3*T*6/2700);                          % From slides, causes high p_plenum (15%)
    %bottom ones are best
    
    AF = mair/mfuel;
    lambda = AF/AFstoi;
    
    EGR = 20;
    EGRf = EGR/100;
    Vd   = pi*(Bore/2)^2*Stroke;                              % Displacement volume for all cylinders
    rho = (1+EGRf)*mair/Vd;
end
% Real AF ratio

fracfu  = 1/(AF+1);
fracair = 1 - fracfu;
YReactants = fracfu*Yfuel + fracair*Yair;                       % Composition vector
YProducts  = YReactants - si*YReactants(1);                     % Corresponding fully combusted setting.
for i=1:Nsp
    hi(i) = HNasa(T0,SpS(i));
end
QLHV = (hi*YReactants'-hi*YProducts')/YReactants(1);            % Classical definition of lower heating value. Just for reference not used!
Comb.QLHV    = QLHV;

Int.T   = T_plenum;
Exh.T   = T_exhaust;
Int.Y   = (1-EGRf)*YReactants+EGRf*YProducts;                   % Applying EGR setpoint
Exh.Y   = YProducts;
Int.T   = T_plenum;
Exh.T   = T_exhaust;
for ii=1:Nsp
    Int.h(ii)   = HNasa(Int.T,SpS(ii));
    Int.e(ii)   = ENasa(Int.T,SpS(ii));
    Int.Cp(ii)  = CpNasa(Int.T,SpS(ii));
    Exh.h(ii)   = HNasa(Exh.T,SpS(ii));
    Exh.e(ii)   = ENasa(Exh.T,SpS(ii));
    Exh.Cp(ii)  = CpNasa(Exh.T,SpS(ii));
end
Mave    = 1/sum([Int.Y]./Mi);
Rg      = Runiv/Mave;

if strcmp(mode,'couple')
    p_plenum = rho*Rg*T_plenum;
%     p_plenum = 3.63e5;
    p_exhaust = p_plenum+0.1*bara;  % exhaust back-pressure
end

% p0      = 3.5*bara;                                             % Typical full load set point
p0 = p_plenum;

mass    = p0*V0/Rg/T0;
massfu  = Int.Y(1)*mass;
Settings.N      = N;
Settings.EGR    = EGRf;
Settings.AF     = AF;
Settings.Ncyc   = Ncyc;             %Added info about the number of cycles



Int.Ca=CaI;Int.L=LI;Int.D=Di;Int.p=p_plenum;
Exh.Ca=CaE;Exh.L=LE;Exh.D=De;Exh.p=p_exhaust;


%% Set initial solution (it is an DAE problem so we must initialize)
y0(1)=p0;y0(2)=T0;y0(3:3+Nsp-1) = mass*[Int.Y];
yNames={'p','T','','','','',''};
for i=3:3+length(Names)-1
    yNames{i}=char(Names(i-2));
end
mfuIVCClose     = y0(3);

%% Save current case to pass on to FtyDAE
if strcmp(mode,'case')
    currentCase = allCases(iCase-121);
    
elseif strcmp(mode,'couple')
    currentCase.T = T;
    currentCase.EGRf = EGRf;
    currentCase.p_plenum = p_plenum;
    currentCase.T_plenum=T_plenum;
    currentCase.w = w;
end

currentCase.mode = mode;
save('currentCase.mat','currentCase');


%% Computing CA values

ReducedCA = -360:360; 
HRR = EtaComb*QLHV*mfuIVCClose*wiebefunctions(ReducedCA,TSOI,PSOI);
% HRR = EtaComb*QLHV*mfuel*wiebefunctions(ReducedCA);

for i = 1:length(HRR)
    %HR(i) = trapz(HRR(1:ReducedCA(i)));
    HR(i) = trapz(HRR(1:i));
end

currentCase.HR = HR;
save('currentCase.mat','currentCase');

% Useful and important Crank Angles determined from heat release rate:
CA01 = ReducedCA(find(HR>0.01*HR(length(HR)),1)+1);
CA05 = ReducedCA(find(HR>0.05*HR(length(HR)),1)+1);
CA10 = ReducedCA(find(HR>0.1*HR(length(HR)),1)+1);
CA50 = ReducedCA(find(HR>0.5*HR(length(HR)),1)+1);
CA90 = ReducedCA(find(HR>0.9*HR(length(HR)),1)+1);
CA95 = ReducedCA(find(HR>0.95*HR(length(HR)),1)+1);
BDUR = CA90-CA01;


%% Solving the DAE system
tspan=t;
odopt=odeset('RelTol',1e-4,'Mass',@MassDAE,'MassSingular','yes');           % Set solver settings (it is a DAE so ...,'MassSingular','yes')
tic;
Qheatloss = 0;      % Total amount of heat lost during the simulation
alfaplot = [];
Vr = max(CylVolumeFie(t));
[time,y]=ode15s(@FtyDAE,tspan,y0,odopt);                                    % Take a specific solver
% disp(Qheatloss);    % Remove comment to display total amount of heat lost for a single complete case
% plot(0:max(time)/(length(alfaplot)-1):max(time),alfaplot); % Remove comment to display alfa values over time
tel=toc;
fprintf('Spent time %9.2f (solver %s)\n',tel,'ode15s');


%% Specify SaveName
if strcmp(mode,'case')
    CaseName = ['Case' num2str(iCase,'%3.3i') '.mat'];
    
elseif strcmp(mode,'couple')
    CaseName = ['paramCase' num2str(paramNumber,'%3.3i') '.mat'];
    
end

% Saves the complete data to the chosen case name:
SaveName = fullfile(DataDir,CaseName);
V = CylVolumeFie(time);
save(SaveName,'Settings','Cyl','Int','Exh','Comb','time','y','yNames','V','SpS','HRR','CA10','CA50','CA90','BDUR','mode','HRR');
fprintf('Saved solution of Case %3i to %s\n',iCase,SaveName);
% save('CaPlots.mat','reducedCaPlot','CaPlot');
end
