function ds = SCM_ODEsolvefunc(t,s,k,Nav,cavpressureparam1, cavpressureparam2, cavpressureparam3, rholiq, nuliq, cwliq, sigmaliq, Dij, lambda, theta, fi, Pv, Tinf, VdW_a_i, VdW_b_i, ncomp, m, eta, paramincludekm, paramincludediff, paramincludetemp, paramincludecvcorr, paramincludeucorr,paramincludediffw,paramincludepsicorr,cavparam,paramterminate)

%% Number of components in the bubble & number of other ODE's being solved 
lengthinit = length(s); % No of initial conditions for the ODE
nodevar = (lengthinit - ncomp); % No of ode's being solved other than that of diffusion
%% 
    
%% Variable for computing the rate change of R, dR/dt, d(components)/dt, dT/dt 
ds = zeros(lengthinit,1);
%% 

%% Assigning initial values into recognizable variable names 
R = s(1);
Rt = s(2);
Tb = s(3);

Ntbubble = sum(s((nodevar+1):lengthinit));

Nbubble = zeros(ncomp,1);
xbubble = zeros(ncomp,1);
VdW_a = 0;
VdW_b = 0;

for i = 1:ncomp
    Nbubble(i) = s(nodevar+i); %Assigning s(4) to s(6) to components in the bubble
    xbubble(i) = (Nbubble(i)/Ntbubble);
    VdW_a = VdW_a + VdW_a_i(i)*(xbubble(i)^2);
    VdW_b = VdW_b + VdW_b_i(i)*xbubble(i);
end

% Ntbubble = sum(Nbubble);

%% Computing values of Rsquare, V, dV/dt for subsequent use in the program 
Rsquare = R*R;
Rtsquare = Rt*Rt;
Rcube = R*R*R;
Rbypi = R/pi;
Vb = ((4*pi/3)*(Rcube));
Vt = (4*pi*Rsquare*Rt);
nbubble = Ntbubble/Nav; % Number of kmol
%%

%% Pressure inside the bubble & bulk medium
% Pb = ((Ntbubble*k*Tb)/(Vb - nbubble*VdW_b));
Pb = ((Ntbubble*k*Tb)/(Vb - nbubble*VdW_b)) - VdW_a*(nbubble^2)/(Vb^2);
Pinf = SCM_Cavpressurefunc(cavparam,cavpressureparam1,cavpressureparam2,cavpressureparam3,t);
%%

%% Bubble calculations 
% xbubble = zeros(ncomp,1);
% for i = 1: ncomp
%     xbubble(i) = (Nbubble(i)/Ntbubble);
% end

Cbubble = zeros(ncomp,1);
for i = 1:ncomp %Concentration of species in the bubble - (Units - no of molecules/m3)
    Cbubble (i) = ((Nbubble(i))/(Vb));
end
%% 

%% Interface calculations 
Pinterface = zeros(ncomp,1); 
Pinterface (1) = (((Pb-Pv)*xbubble(1)));%Partial pressure of species at the interface - (Units - Pa)
Pinterface (2) = (((Pb-Pv)*xbubble(2)));%Partial pressure of species at the interface - (Units - Pa)
Pinterface (3) = Pv; %Partial pressure of species at the interface - (Units - Pa)

Cinterface = zeros(ncomp,1);
for i = 1: ncomp %Concentration of species at the interface- (Units - no of molecules/m3)
    Cinterface (i) = ((Pinterface(i))/ (k*Tinf));
end

xinterface = zeros(ncomp,1);
Ctinterface = sum(Cinterface);

for i = 1:ncomp %mole fraction of species at interface - (Units - (-))

    xinterface(i) = ((Cinterface(i))/(Ctinterface));
end

Cinterfaceij = zeros(3,3);

for i = 1: ncomp 
    for j = 1: ncomp 
        Cinterfaceij(i,j) = (Cinterface (i) + Cinterface (j));
    end
end
%% 

%% 

%% Diffusivity calculation 

dbynij = zeros(ncomp,ncomp); 
for i = 1 : ncomp  
    for j = 1 : ncomp 
        dbynij (i,j) = (Dij(i,j)/(Cinterfaceij(i,j)));
    end
end

Di = zeros(ncomp,1);
Diinverse = zeros(ncomp,1);

for i = 1 : ncomp % Di inverse calculation - (1/Diffusivity) (Units -s2/m)
    for j = 1: ncomp
        if (i~=j)
            Diinverse(i) = (Diinverse(i) + ((xbubble(j))/ ((1-xbubble(i)) * (dbynij(i,j)))));
        end
    end
end 

for i = 1: ncomp 
    Di(i) = (1/Diinverse(i));
end 
%% 

%% Thermal conductivity of contents in the bubble
lmbdmx = 0; 
phi(ncomp,ncomp) = 0;
for i = 1:ncomp 
    
    denom = 0;
      
    for j = 1:ncomp 
        
        phidummy1 = m(i);
        phidummy2 = m(j);
        phidummy3 = (phidummy1/phidummy2);
                
        phidummy4 = eta(i);
        phidummy5 = eta(j);
        phidummy6 = (phidummy4/phidummy5);
        
        phidummy7 = (8*(1+phidummy3));
        phidummy8 = (sqrt(phidummy7));
        
        phidummy9 = ((phidummy6)^-0.5); %eta term
        phidummy10 = ((phidummy3)^0.25); %m term
        phidummy11 = (1 +(phidummy9*phidummy10))^2;
        
        phi(i,j) = phidummy11/phidummy8;
        denom = (denom + (xbubble(j)*phi(i,j)));                
        
               
    end
    
    lmbdmx = (lmbdmx + ((xbubble(i)*lambda(i))/(denom)));
    
end 
%% 

%% Cv calculation
Cv = zeros(ncomp,1);
Psi = zeros(ncomp,1);
for i = 1: ncomp
    Cv(i) = SCM_Cvcalc(k,fi(i),theta(i,:),Tb,paramincludecvcorr);
    Psi(i) = SCM_Psicalc(k,theta(i,:),Tb,paramincludepsicorr);
end

% for i = 1: ncomp
%     for j = 1:length(theta)
%         Cv(i) = SCM_Cvcalc(k,fi(i),theta(i,j),Tb,paramincludecvcorr);
%     end
% end

Cvmix = 0;

for i = 1: ncomp  
    Cvmix = (Cvmix + (Nbubble(i)*(Cv(i) - Psi(i))));
end


%% 
    

%% Cp & kappa calculation
Cp = zeros(ncomp,1);
for i = 1: ncomp
    temp1 = (fi(i) + 2);
    temp2 = (temp1/2);
    temp3 = (temp2 * k);
    Cp(i) = (Cp(i) + temp3);
end

rhocpmix = zeros(ncomp,1);
for i = 1:ncomp
    rhocpmix (i) = ((Cp(i))*(Cbubble(i)));
end

rhocpmixsum = sum(rhocpmix);

kappa = (lmbdmx/ (rhocpmixsum));
%% 

%% ldiff & lth calculation
ldiff = zeros(ncomp,1); %calculation of ldiff - (Units - m2/s)

for i = 1: ncomp
    ldiff(i) = sqrt ((R*Di(i))/abs(Rt));
end
    
lth = sqrt ((R*kappa)/abs(Rt));

for i = 1:ncomp
    if (ldiff(i) > Rbypi)
        ldiff(i) = Rbypi;
    end
end


if(lth > Rbypi)
    lth = Rbypi;
end
%% 


%% Rate change of species calculation
for i = (nodevar+1) : (lengthinit-1) %Rate change for species other than water...
     ds (i) = ((4*pi*Rsquare)* (Di(i-3))* paramincludediff * ((Cinterface(i-3)-Cbubble(i-3))/ldiff(i-3)));
end


%Rate change for water...
ds(lengthinit) = ((4*pi*Rsquare)* (Di(lengthinit-nodevar))* paramincludediffw * ((Cinterface(lengthinit-nodevar)-Cbubble(lengthinit-nodevar))/ldiff(lengthinit-nodevar)));

dNbydt = 0;

 for i = (nodevar+1):(lengthinit)
    dNbydt = dNbydt + ds(i);
 end
 %% 

%% Enthalphy & Internal energy calculation
hi = zeros(ncomp,1);
for i = 1:ncomp
    hi (i) = (Cp(i)*Tinf);
end

Ui = zeros(ncomp,1);

for i = 1:ncomp
    Ui (i) = SCM_Uifunc(fi(i),k,Tb,theta(i,:),paramincludeucorr);
end

%% Components in the energy balance equation
dNprod = zeros(ncomp,1);

for i = 1: ncomp
%     dNprod (i) = ((hi (i) - Ui (i))* ds(i+3));
%     dNprod (i) = ((hi (i) - Ui (i) - Tb*Cv(i))* ds(nodevar+i));
    dNprod (i) = ((hi (i) - Tb*Cv(i))* ds(nodevar+i));

end

dNprodsum = sum(dNprod);

Qtterm1 = (4*pi*Rsquare);
Qtterm2 = (Tinf - Tb);
Qtterm3 = (lmbdmx);
Qt = ((Qtterm1)*(Qtterm2/lth)*(Qtterm3)); %Splitting the terms of Qt; 3rd Feb 2018 - Varaha

Pbvt = (Pb*Vt);

if (paramincludetemp == 1)
    dTbydt = ((Qt - Pbvt + dNprodsum)/ Cvmix);
else
    dTbydt = 0;
end
%% 


%% Inclusion/Exclusion of correction in Keller-Miksis 
if (paramincludekm == 1)
    phi1 = (1 - (Rt/cwliq));
else
    phi1 = 1;
end
    
if (paramincludekm == 1)
    phi2 = (1 - (Rt/(3*cwliq)));
else
    phi2 = 1;
end

if (paramincludekm == 1)
    phi3 = (1 + (Rt/cwliq));
else
    phi3 = 1;
end
%% 
        
%% Terms for the calculation of (dpb/dt) - Version March 2020
% dpbydtterm1 = ((Tb*dNbydt)+(Ntbubble*dTbydt));
% dpbydtterm2 = (3*Rsquare*Rt*Ntbubble*Tb);
% dpbydtterm3 = (Rcube - hccube);
% dpbydtterm4 = ((dpbydtterm1)/(dpbydtterm3));
% dpbydtterm5 =  (dpbydtterm3*dpbydtterm3);
% dpbydtterm6 = ((dpbydtterm2)/(dpbydtterm5 ));
% dpbydtterm7 = (dpbydtterm4 - dpbydtterm6);
% dpbydt = ((3*k)/(4*pi))*dpbydtterm7;
%% Terms for the calculation of (dpb/dt) - Version June 2020
% dpbydtterm1 = (Vb - nbubble*VdW_b)*(k*Tb*dNbydt + k*Ntbubble*dTbydt);
% dpbydtterm2 = Ntbubble*k*Tb*(Vt - (VdW_b/Nav)*dNbydt);
% dpbydtterm3 = (Vb - nbubble*VdW_b)^2; 
% dpbydt = (dpbydtterm1 - dpbydtterm2)/dpbydtterm3;
%% Terms for the calculation of (dpb/dt) - Version August 2020
dpbydtterm1 = (Vb - nbubble*VdW_b)*(k*Tb*dNbydt + k*Ntbubble*dTbydt);
dpbydtterm2a = Ntbubble*k*Tb*Vt;
dpbydtterm2b = 0;
dpbydtterm3 = (Vb - nbubble*VdW_b)^2; 
dpbydtterm4 = (6*VdW_a*(nbubble^2)*Rt)/(R*(Vb^2));
dpbydtterm5 = 0;
for i = 1:ncomp
    dpbydtterm2b = dpbydtterm2b + VdW_b_i(i)*ds(nodevar+i)/Nav;
    dpbydtterm5 = dpbydtterm5 + (2*VdW_a_i(i)/Nav)*((Nbubble(i)/Nav)/(Vb^2))*ds(nodevar+i);
end
dpbydtterm2b = Ntbubble*k*Tb*dpbydtterm2b;
dpbydtterm2 = dpbydtterm2a - dpbydtterm2b;


dpbydt = (dpbydtterm1 - dpbydtterm2)/dpbydtterm3 + dpbydtterm4 - dpbydtterm5;
%% 

%% Terms for the calculation of (dR/dt)
Rtterm1 = ((phi3*(Pb-Pinf))/(phi1*rholiq*R));
Rtterm2 = ((dpbydt)/(phi1*rholiq*cwliq));
Rtterm3 = ((4*nuliq*Rt)/(phi1*Rsquare));
Rtterm4 = ((2*sigmaliq)/(phi1*rholiq*Rsquare));
Rtterm5 = ((1.5*phi2*Rtsquare)/(phi1*R));
%%

%% Assigning calculated values to the ODE
ds(1) = Rt;
ds(2) = (Rtterm1 + Rtterm2 - Rtterm3 - Rtterm4 - Rtterm5);
ds(3) = dTbydt;
%% 

end