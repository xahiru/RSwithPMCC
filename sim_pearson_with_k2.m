function [ sims ] = sim_pearson_with_k2( prefs,k )
%sim_pearson_with_k2( prefs,k ) calculates pairwise correlations between
%the raws of prefs matrix with threshold of k
%   Detailed explanation goes here

% sims = zeros(size(prefs,1),size(prefs,1));

  sims = nan((size(prefs,1)),(size(prefs,1)));
%  sims = nan(7,7);

for j=1:size(prefs,1)
        
      for i=1:size(prefs,1)
        
%              if i~=j

%         if(k == 6)
%             sims(j,i)  = sim_own(prefs,i,j,k);
%         else
            
            sims(i,j)  =  sim_pearson_with_k(prefs,i,j,k);
     
%              end
        
%    sims(j,i)  = sim_pearson(prefs,i,j,k);
         
           
%            end
    
     end



end
%      sims(logical(eye(size(sims)))) = [];
%       sims = reshape(sims,size(prefs,1),size(prefs,1)-1);
%     k = 

end
