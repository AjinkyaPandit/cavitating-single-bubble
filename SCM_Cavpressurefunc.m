%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to estimate the pressure profile of Acoustic/Hydrodynamic Cavitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pt = SCM_Cavpressurefunc(cavparam,cavpressureparam1,cavpressureparam2,cavpressureparam3,t)

   
if(cavparam == 1) %Pt = [P0 - (Pa * sin (omega*t))]
    P0 = cavpressureparam1;
    Pa = cavpressureparam2;
    omega = cavpressureparam3;
    f = omega/(2*pi);
%     P0 = t*P0;
    Pt = (P0 - (Pa * sin((omega*t))));
%     if (f*t > 3)
%         Pt = P0;
%     end
         
%% 

%APCRE input - October 2017

elseif(cavparam == 2)

    PTh = cavpressureparam1;
    P2 = cavpressureparam2;
    Q = cavpressureparam3;

    Tau = 1.6e-3;
    rhoL = 1000;
    uTh = 29.12; 
    uFl = 2.55;
    omega = 4.55*1000*2*pi;

    Pt = PTh + ((P2 - PTh)*t)/(Tau);
    ut = sqrt((PTh + 0.5*rhoL*(uTh^2) - Pt)/(0.5*rhoL));
    ut_new = ut + uFl*sin(omega*t);
    Pt = PTh + 0.5*rhoL*(uTh^2) - 0.5*rhoL*(ut_new^2);

elseif(cavparam == 3)
    td = cavpressureparam1;
    Ptd = cavpressureparam2;
    
    %Pt = spline(cavpressureparam1,cavpressureparam2,t);
    Pt = pchip(cavpressureparam1,cavpressureparam2,t);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%