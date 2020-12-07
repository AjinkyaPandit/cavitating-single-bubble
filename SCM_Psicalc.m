function Psii = SCM_Psicalc(k,theta,Tb,paramincludepsicorr)
    
    Psiisum = 0;
    Psiisum1 = 0;
    
%         if(theta~=0)
%             var = (theta/(Tb));
%             tempsq = var*var;
%             Cvsum1 = (((tempsq)*exp(var))/((exp(var)-1)^2));
%             
%         end

   for j = 1: length (theta)
        if(theta(j)~=0)
            var = (theta(j)/(Tb));
            varsq = var*var;
            psiterm1 = (exp(var)^2)*(2*var - varsq);
            psiterm2 = (varsq + 2*var)*(exp(var));
            psiterm3 = ((exp(var) - 1)^3);
            Psiisum1 = Psiisum1 + var*(psiterm1 - psiterm2)/(psiterm3);
        end
    end
                
   if (paramincludepsicorr == 1)
       Psiisum = Psiisum1;
   else
       Psiisum = 0;
   end
    
%    Cv = (Cvsum*k*Tb);
   Psii = (Psiisum*k);
   
end
