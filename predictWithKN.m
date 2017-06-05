function [ Predicted ] = predictWithKN( user,prefs,KN, simsk2u, itemLength )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


kevin_corrs = simsk2u(user,:);
kevin_corrs(user) = [];

[a, b] = sort(kevin_corrs,'descend');


 kevin_mu = nanmean(prefs(user,:));

ngh_sim = a(1:KN);
ngh_idx = b(1:KN);

 
%  ngh_sim(invKN:end) = [];                       % drop non-neighbors
%  ngh_idx(invKN:end) = [];                        % drop non-neighbors
% ngh_sim((length(ngh_sim)-(KN-1)):end) = [];
% ngh_idx((length(ngh_idx)-(KN-1)):end) = [];

ngh_mu = nanmean(prefs(ngh_idx,:));
Predicted = nan(itemLength,1);

for i = 1:itemLength                    % loop over movies
    ngh_r = prefs(ngh_idx,i);             % neighbor ratings for the movie
%      isRated = any(ngh_r);                % only use neighbors who rated
       isRated = ~isnan(ngh_r); 
    meanCentered =...                       % mean centered weighted average
        (ngh_r(isRated)' - ngh_mu(isRated)) * ngh_sim(isRated)'...
        / sum(ngh_sim(isRated));
 

 	if isempty(meanCentered)
     meanCentered = NaN;
    Predicted(i) =  NaN;
    else
        Predicted(i) =  kevin_mu + meanCentered; % add Kevin's average
     end

    
end





end

