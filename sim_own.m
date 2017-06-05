function [ sim ] = sim_own( prefs,x1,x2,k )
%Simlarity between coloumn vector x1 and x2 from preference matrix
%peasron correlation sim value is returned
%   Detailed explanation goes here

    c1=find(prefs(x1,:)>0);
    c2=find(prefs(x2,:)>0);
    
    cc=intersect(c1,c2);
    
    shared_items(1,:)=prefs(x1,cc);
    shared_items(2,:)=prefs(x2,cc);
    
   
    
    if size(shared_items,2)== 0
        sim = 0;
    else
        % if not, run the algorithm
        x=shared_items(1,:);
        y=shared_items(2,:);
        
        
        xmean = mean(x);
        ymean = mean(y);
        
        num = sum((x-xmean).*(y-ymean));
        
        
         den = sqrt(sum((x-xmean).^2))*sqrt(sum((y).^2));
         
         if den==0
            sim = 0;
        else
         
         sim = num/den;
         end
    end

end

