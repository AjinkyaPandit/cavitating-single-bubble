%% Single Cavity Main Program 

run SCM_ReadData
run SCM_estimatederivedprop

% init = [(2^(1/3))*7.70E-05, 0.0, 301.611132, 2*25592213508, 2*96275469863, 2*1.58934E+12];
% time = [7.40E-05, time(2)];
[t,s] = ode15s(@SCM_ODEsolvefunc, time, init, options, k, Nav, cavpressureparam1, cavpressureparam2, cavpressureparam3, rholiq, nuliq, cwliq, sigmaliq, Dij, lambda, theta, fi, Pv, Tinf, VdW_a_i, VdW_b_i, ncomp, m, eta, paramincludekm, paramincludediff, paramincludetemp, paramincludecvcorr, paramincludeucorr,paramincludediffw,paramincludepsicorr,cavparam,paramterminate);

run SCM_Postprocess
