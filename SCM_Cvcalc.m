function Cv = SCM_Cvcalc(k,fi,theta,Tb,paramincludecvcorr)
    
    Cvsum = 0;
    Cvsum1 = 0;
    
%         if(theta~=0)
%             var = (theta/(Tb));
%             tempsq = var*var;
%             Cvsum1 = (((tempsq)*exp(var))/((exp(var)-1)^2));
%             
%         end

   for i = 1: length (theta)
        if(theta(i)~=0)
            var = (theta(i)/(Tb));
            tempsq = var*var;
            Cvsum1 = Cvsum1 + (((tempsq)*exp(var))/((exp(var)-1)^2));
        end
    end
                
   if (paramincludecvcorr == 1)
       Cvsum = (fi/2) + Cvsum1;
   else
       Cvsum = (fi/2);
   end
    
%    Cv = (Cvsum*k*Tb);
   Cv = (Cvsum*k);
   
end
