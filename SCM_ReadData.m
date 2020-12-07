%% SCM_Constants - Avagadro & Boltzmann 
Nav = 6.02e26; % Avagadro constant, (no. of molecules/gmol)
k = 1.38e-23; % Boltzmann constant, (J/(K * no of molecules)
%% 

%% SCM_Dissolvedgasproperties
fi = xlsread('InputforSCM.xlsx','Dissolvedgasprop','B2:B4'); 
Mw = xlsread('InputforSCM.xlsx','Dissolvedgasprop','C2:C4'); 
alpha = xlsread('InputforSCM.xlsx','Dissolvedgasprop','D2:D4'); 
epsibyk= xlsread('InputforSCM.xlsx','Dissolvedgasprop','E2:E4');
theta = xlsread('InputforSCM.xlsx','Dissolvedgasprop','F2:L4'); 
ncomp = length(fi);
%% 

%% SCM_Liquidproperties
rholiq = xlsread('InputforSCM.xlsx','Liquidprop','B2');
muliq = xlsread('InputforSCM.xlsx','Liquidprop','B3');
sigmaliq = xlsread('InputforSCM.xlsx','Liquidprop','B4');
cwliq = xlsread('InputforSCM.xlsx','Liquidprop','B5');
Pv = xlsread('InputforSCM.xlsx','Pressureprofile','B19');

%% SCM_Cavitation Device & Pressure
cavparam = xlsread('InputforSCM.xlsx','Pressureprofile','B14');

if(cavparam == 1)
     cavpressureparam1 = xlsread('InputforSCM.xlsx','Pressureprofile','B5');
     cavpressureparam2 = xlsread('InputforSCM.xlsx','Pressureprofile','B6');
     freq = xlsread('InputforSCM.xlsx','Pressureprofile','B7');
     freq = freq*1000;
     cavpressureparam3 = 2*pi*freq;
      
else if(cavparam == 2)
     cavpressureparam1 = xlsread('InputforSCM.xlsx','Pressureprofile','B10');
     cavpressureparam2 = xlsread('InputforSCM.xlsx','Pressureprofile','B11');
     cavpressureparam3 = xlsread('InputforSCM.xlsx','Pressureprofile','B12');
     
    else
        
     cavpressureparam1 = xlsread('InputforPressure.xlsx','PressureValues', 'A2:A5002');
     cavpressureparam2 = xlsread('InputforPressure.xlsx','PressureValues', 'B2:B5002');
     cavpressureparam3 = 1; %dummy variable in case of interpolating pressure
    end
end 
%% 

%% SCM_Correctionfactors

paramincludekm = xlsread('InputforSCM.xlsx','Modelparam','B20');
paramincludediff = xlsread('InputforSCM.xlsx','Modelparam','B21');
paramincludetemp = xlsread('InputforSCM.xlsx','Modelparam','B22');
paramincludecvcorr = xlsread('InputforSCM.xlsx','Modelparam','B23');
paramincludeucorr = xlsread('InputforSCM.xlsx','Modelparam','B24');
paramincludediffw = xlsread('InputforSCM.xlsx','Modelparam','B25');
paramincludepsicorr = xlsread('InputforSCM.xlsx','Modelparam','B26');
paramterminate = xlsread('InputforSCM.xlsx','Modelparam','B27');
%% 

%% Ambient pressure, temperature and initial conditions for the simulation
Tinf = xlsread('InputforSCM.xlsx','Pressureprofile','B1');
Pinf = xlsread('InputforSCM.xlsx','Pressureprofile','B5');
init = xlsread('InputforSCM.xlsx','Modelparam','B2:B9');
%% 

%% Initial radius of bubble and Van der Waal's hard core radius
R0 = init(1);
VdW_a_i = xlsread('InputforSCM.xlsx','Modelparam','I30:I32');
VdW_b_i = xlsread('InputforSCM.xlsx','Modelparam','B30:B32');

%% Reading time span for the simulation
tinit = xlsread('InputforSCM.xlsx','Modelparam','B10');

if(cavparam == 1)
tend = xlsread('InputforSCM.xlsx','Modelparam','B11');
else
tend = xlsread('InputforSCM.xlsx','Modelparam','B12');
end

time = [tinit tend];
%% 

%% MATLAB solver tolerances
abstol= xlsread('InputforSCM.xlsx','Modelparam','B28');
reltol= xlsread('InputforSCM.xlsx','Modelparam','B29');
% options = odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn', @odeplot,'Events',@SCM_event_function);
options = odeset('RelTol',reltol,'AbsTol',abstol,'Events',@SCM_event_function);