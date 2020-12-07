%% Vapor pressure & kinematic viscosity of water
nuliq = muliq/rholiq; %Kinematic viscosity of water, (m2/s)
%% 

%% Number of components in the bubble & number of other ODE's being solved 
lengthinit = length(init); % No of initial conditions for the ODE
nodevar = (lengthinit - ncomp); % No of ode's being solved other than that of diffusion
%%

%% (Dij * nij ) & lambda calculation
Tred1 = zeros(ncomp,1);
omega22 = zeros(ncomp,1);
m = zeros(ncomp,1);
eta = zeros(ncomp,1);
lambda = zeros(ncomp,1);

Tred2 = zeros(ncomp,1);
omega11 = zeros(ncomp,1);
tau = zeros(ncomp,1);
deff = zeros(ncomp,1);
Dij = zeros(ncomp,ncomp);
VdW_a = 0;
VdW_b = 0;

for i = 1:ncomp
    
    Tred1 (i) = (Tinf/epsibyk(i));
    omega22 (i) = (-34.90317 + (0.83841*(Tred1(i)^-0.73653)) + (35.6314*(Tred1(i)^-0.00109)));
    m(i) = (Mw(i)/(Nav));
    eta(i) = (5/16)*sqrt(pi*m(i)*k*Tinf)/(pi*alpha(i)^2*omega22(i));
    lambda(i) = (15/4)*(k/m(i))*eta(i)*((4/15)*(fi(i)/2) + (3/5)); %Thermal conductivity
    VdW_a = VdW_a + VdW_a_i(i)*(init(nodevar+i)^2);
    VdW_b = VdW_b + VdW_b_i(i)*init(nodevar+i);
    
    for j = 1:i
        
        Tred2(i,j) = (Tinf/(sqrt(epsibyk(i)*epsibyk(j))));
        omega11(i,j) = [(-34.8998 + (0.69463*(Tred2(i,j)^-0.825232)) + (35.63592*(Tred2(i,j)^-0.00146)))];
        tau(i,j) = ((2*m(i)*m(j))/(m(i) + m(j)));
        deff(i,j) = 0.5*(alpha(i) + alpha(j));     
        Dij(i,j) = ((3/8)*sqrt(pi*k*Tinf/tau(i,j)))/((pi*deff(i,j)^2)*omega11(i,j));
        Dij(j,i) = Dij(i,j);
        
    end
    
end
%% 


%% Initial - pressure, volume and number of molecules in the bubble
if(cavparam == 1) % for AC 
    Pbubble = (Pinf + ((2*sigmaliq)/R0));
else % for HC
    Pbubble = (Pv + ((2*sigmaliq)/R0));
end

%% VdW Volume Correction
% Vbubble = (4*pi/3)*(R0*R0*R0);
Vbubble = (4*pi/3)*((R0^3));
Ninit = ((Pbubble*Vbubble)/((Pbubble*VdW_b/Nav) + k*Tinf)); %Number of molecules at time t = 0 (Units - no of molecules)

%% VdW Pressure Correction
% fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
fun = @(x)Pbubble*Vbubble - (VdW_b*Pbubble + k*Nav*Tinf)*x - VdW_a*((x^2)/Vbubble) + VdW_a*VdW_b*((x^3)/(Vbubble^2));
x0 = Ninit/Nav;
x = fminsearch(fun,x0);
Ninit = x*Nav;

%% Ideal Gas
% Ninit = Pbubble*Vbubble/(k*Tinf);
%% 

%% Calculation of initial number of molecules in the bubble from initial mole fractions
% init (ncomp + i) = 0;
% nstate = length(init) - ncomp; % Number of ODEs for state variables
for i = 1:ncomp
    init(nodevar+i) = Ninit*init(nodevar+i);
end