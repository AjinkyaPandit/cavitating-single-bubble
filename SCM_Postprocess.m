%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file has been written to post - process results from the SCM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots in MATLAB
[numb ~] = size(t);
Ntot = (s(:,4) + s(:,5) + s(:,6));
Nwater = s(:,6);
xw = (Nwater./Ntot);

radarray = (s(:,1).^3);

Varray = ((4*pi/3).*radarray);
Tarray = s(:,3);

xBubble = zeros(numb,ncomp);
VdW_a = zeros(numb,1);
VdW_b = zeros(numb,1);
for j = 1:numb
    for i = 1:ncomp
        xBubble(j,i) = s(j,nodevar+i)/Ntot(j);
        VdW_a(j) = VdW_a(j) + VdW_a_i(i)*(xBubble(j,i)^2);
        VdW_b(j) = VdW_b(j) + VdW_b_i(i)*xBubble(j,i);
    end
end

% a = 5/101325;
% Pmaxarray = ((Ntot.*Tarray*k)./(Varray - Ntot.*(VdW_b./Nav))) - a/(Varray.*Varray);

Pmaxarray = ((Ntot.*Tarray*k)./(Varray - Ntot.*(VdW_b./Nav))) - VdW_a.*(Ntot.*Ntot/(Nav^2))./(Varray.*Varray);
Pmaxarray = Pmaxarray/101325;

Pt = zeros(numb,1);
for i=1:numb
    Pt(i) = SCM_Cavpressurefunc(cavparam,cavpressureparam1,cavpressureparam2,cavpressureparam3,t(i));
end
subplot(2,2,1)       
plot(t,s(:,1)/R0)           
xlabel('Time, (s)');
ylabel('R/R0, (-)');
title('Relative Radius of Bubble');
% yyaxis right;
% plot(t,Pt);

subplot(2,2,2)  
plot(t,s(:,3))
xlabel('Time, (s)');
ylabel('T, (K)');
title('Temperature inside the bubble');

subplot(2,2,3)
plot(t,Pmaxarray, 'b', t, Pt/101325, 'r');
xlabel('Time, (s)');
ylabel('Pressure, (bar)');
title('Pressure inside the bubble');

subplot(2,2,4)  
plot(t,xw)
xlabel('Time, (s)');
ylabel('Mole fraction of water');
title('xw');

[Tmax,loc]= max(Tarray);

% output = [t(loc) max(Tarray) max(Pmaxarray) xw(loc)];
% output=output';


%% 

% if (cavparam == 1) %Pt = [P0 - (Pa * sin (omega*t))]
%     P0 = cavpressureparam1;
%     Pa = cavpressureparam2;
%     omega = cavpressureparam3;
%     Pt = (P0 - (Pa * sin((omega*t))));
%          
% else
%      td = cavpressureparam1;
%      Ptd = cavpressureparam2;
%      
%      Pt = pchip(cavpressureparam1,cavpressureparam2,t);
%      
% end
% 
%% Bubble Energy Calculations

Eb = zeros(numb,1);
Cv = zeros(ncomp,1);
for i = 1:numb
    % Cv calculation
    Tb = s(i,3);
    Nb = s(i,nodevar+1:nodevar+ncomp);
    for j = 1: ncomp
        Cv(j) = SCM_Cvcalc(k,fi(j),theta(j,:),Tb,paramincludecvcorr);
        Eb(i) = Eb(i) + Cv(j)*Nb(j)*Tb;
    end
end

%% Performance Parameters
tplot = t.*freq;
% output = [tplot,s(:,1)*1e6,Tarray,Pmaxarray,Nwater,Pt];
% output = [max(Pmaxarray), max(Tarray), max(s(:,1)/s(1,1)), max(Nwater)];
% output = [tplot, t, s, Pmaxarray];
% output = [tplot, s, Pmaxarray];
Rmax = max(s(:,1));
v = s(:,2).*(s(:,1).*s(:,1))/(Rmax^2);
v_avg = trapz(t,abs(v))/t(numb);
tke = 0.5*(v.*v);
Iturb = trapz(t,s(:,1))/t(numb);
tke_avg = trapz(t,tke)/t(numb);
eps = (tke.^1.5)/Iturb;
avg_eps = trapz(t,eps)/t(numb);
delP = sqrt(1000*(1e-3)*avg_eps);
beta = cavpressureparam2/((cavpressureparam3^2)*rholiq*(s(1,1)^2));
% output = [delP, max(Pmaxarray), max(s(:,3))];
% output = [delP, max(Pmaxarray), max(s(:,3)), avg_eps, Iturb, k_avg, Rmax];
% [Pmax, ind] = max(Pmaxarray);
% output = [s(1,1)*1e6, Pmax, s(ind,3), s(ind,4), s(ind,5), Nwater(ind)];
% output = [beta, s(1,1)*1e6, delP, max(Pmaxarray), max(s(:,3)), avg_eps, Iturb, tke_avg, Rmax];
nu = cavpressureparam2/cavpressureparam1;
output = [beta, nu, Pmaxarray(numb), s(numb,3), Rmax, Iturb, tke_avg, Rmax];

% output = [t, s, Pmaxarray];

copy(output);

