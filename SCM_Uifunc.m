function Ui = SCM_Uifunc(fi,k,Tb,theta,paramincludeucorr)
    
    Uwsum = 0;
    
    for i = 1:length(theta)
        
        if(theta(i)~=0)
        temp = (theta(i)/Tb);
        Uwsum = (Uwsum + ((temp/(exp(temp)-1))));
        end
        
    end
    
    if(paramincludeucorr == 1)
        Ui = (((fi/2) + Uwsum) * (k*Tb));
    else
        Ui = ((fi/2)*k*Tb);
    end
     
end

        
