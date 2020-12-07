%% Single Cavity Main Program 
run Loop_ReadData;
paramterminate = 1; % Set by Default to one 
init = zeros(3+ncomp,1);
% % paramLoop = 4*Pv:Pv:100*Pv;
paramLoop = 0.0:0.01:2.5;
numLoop = length(paramLoop);
% numLoop = 1;
outputLoop = zeros(numLoop, 9);
for lc = 1:numLoop
    clc;
    P_a = paramLoop(lc)
%     cavpressureparam1 = paramLoop(lc)
%     Pinf = cavpressureparam1;
    cavpressureparam2 =  P_a*cavpressureparam1;
    
    init(1) = R0;
    init(2) = initRead(2);
    init(3) = initRead(3);
    init(4) = initRead(4);
    init(5) = initRead(5);
    init(6) = initRead(6);
    run Loop_estimatederivedprop;

    % init = [(2^(1/3))*7.70E-05, 0.0, 301.611132, 2*25592213508, 2*96275469863, 2*1.58934E+12];
    % time = [7.40E-05, time(2)];
%     [t,s] = ode15s(@SCM_ODEsolvefunc, time, init, options, k, Nav, cavpressureparam1, cavpressureparam2, cavpressureparam3, rholiq, nuliq, cwliq, sigmaliq, Dij, lambda, theta, fi, Pv, Tinf, VdW_b, ncomp, m, eta, paramincludekm, paramincludediff, paramincludetemp, paramincludecvcorr, paramincludeucorr,paramincludediffw,cavparam);
    [t,s] = ode15s(@SCM_ODEsolvefunc, time, init, options, k, Nav, cavpressureparam1, cavpressureparam2, cavpressureparam3, rholiq, nuliq, cwliq, sigmaliq, Dij, lambda, theta, fi, Pv, Tinf, VdW_a_i, VdW_b_i, ncomp, m, eta, paramincludekm, paramincludediff, paramincludetemp, paramincludecvcorr, paramincludeucorr,paramincludediffw,paramincludepsicorr,cavparam,paramterminate);

    run Loop_Postprocess;
    output = [beta, P_a, Pmaxarray(numb), s(numb,3), Rmax, Eb(numb), Rmax, delW_avg, viscE_avg];
    outputLoop(lc,:) = output;
end
  
plot(outputLoop(:,2),outputLoop(:,3))           
xlabel('Time, (s)');
ylabel('R/R0, (-)');
title('Maximum Pressure');

copy(outputLoop);

beep