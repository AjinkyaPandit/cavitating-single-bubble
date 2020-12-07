function [value,isterminal,direction] = SCM_ODEsolvefunc(t,s,k, Nav, cavpressureparam1, cavpressureparam2, cavpressureparam3, rholiq, nuliq, cwliq, sigmaliq, Dij, lambda, theta, fi, Pv, Tinf, VdW_a_i, VdW_b_i, ncomp, m, eta, paramincludekm, paramincludediff, paramincludetemp, paramincludecvcorr, paramincludeucorr,paramincludediffw,paramincludepsicorr,cavparam, paramterminate)
    
value = s(2);  % when value = 0, an event is triggered
% To suppress detection of initial perturbation
f = cavpressureparam3/(2*pi);
tcutoff = 0.125/f;
if t > tcutoff
    isterminal = paramterminate; % terminate after the first event
else
    isterminal = 0;
end
% isterminal = 0;
direction = 1;  % get all the zeros

end