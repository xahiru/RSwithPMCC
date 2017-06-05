%% 
% matlab online CF tutorial
%% load csv into matrix
% filename = '/Users/xahiru/Documents/MATLAB/ml-latest-small/ratings.csv';
filename = '/Users/xahiru/Documents/MATLAB/ml-latest-small/ratings2.csv';
DATA = csvread(filename,1,0);
save('DATA');

clear filename;
% sort data
% C is the user
% B is items
[USER,ia,ic] = unique(DATA(:,1),'rows');
[ITEM,ik,ib] = unique(DATA(:,2),'rows');

%ia is the index of unique id users, ic is the index of all ids the users
%ik is the index of unique id moives, ic is the index of all ids the movies

% create the zero array to hold data
%  mldata=NaN(size(ia,1),size(ik,1));
 mldata=zeros(size(ia,1),size(ik,1));
% go through all the rows to create the data matrix
for j=1:size(DATA(:,1))
    user= find(USER(:,1)== DATA(j,1));
    movie= find(ITEM(:,1)== DATA(j,2));
    rating= DATA(j,3);
    mldata(user,movie)=rating;
    
end

clear rating movie user i

       
 %%
 figure
subplot(2,1,1)
a = 1;
b = 4;
c = 2;
d = 7;

scatter(mldata(a,:),mldata(b,:),'filled')
lsline
xlim([0 5]); ylim([0 5])
title('Movie Preference Space by Two Users')
xlabel(sprintf('User %d', a)); ylabel(sprintf('User %d', b))

for i = 1:size(mldata,2)
%     mldata(i,1)+0.05,mldata(i,2),ITEM(i)
   
    text(mldata(a,i)+0.05,mldata(b,i),int2str(ITEM(i)))
end
subplot(2,1,2)
scatter(mldata(c,:),mldata(d,:),'filled')
lsline
xlim([0 5]); ylim([0 5])
xlabel(sprintf('User %d', c)); ylabel(sprintf('User %d', d))
for i = 1:size(mldata,2)
%     mldata(i,1)+0.05,mldata(i,4),ITEM(i)
    text(mldata(c,i),mldata(d,i),int2str(ITEM(i)))
end


%%
sims = corr(mldata', 'rows', 'pairwise');
% sims = corrcoef(mldata', 'rows', 'pairwise');
fprintf('Similarity between Kevin and Jay:     %.2f\n',sims(2,7))
fprintf('Similarity between Kevin and Spencer:  %.2f\n',sims(2,4))


% normal correlation/similarity user / items
% sims = corr(mldata, 'rows', 'pairwise');
% fprintf('User sim\n')
% fprintf('User Similarity between 1 and 2:  %.2f\n',sims(2,7))
% fprintf('User Similarity between 1 and 4:  %.2f\n',sims(1,4))
% 
% sims = corr(mldata', 'rows', 'pairwise');
% fprintf('Item sim\n')
% 
% fprintf('Item Similarity between 1 and 2:  %.2f\n',sims(2,7))
% fprintf('Item Similarity between 1 and 4:  %.2f\n',sims(1,4))

sims = sims - eye(length(USER)); % set self-correlations to 0
kevin_corrs = sims(4,:);
[ngh_corr, ngh_idx] = sort(kevin_corrs,'descend');
ngh_corr

%%  correlation/similarity with k user / items
% sims = corr(mldata, 'rows', 'pairwise');

simsk2u0 = sim_pearson_with_k2(mldata,0);
simsk2u = sim_pearson_with_k2(mldata,1);
simsk2u2 = sim_pearson_with_k2(mldata,2);
simsk2u3 = sim_pearson_with_k2(mldata,3);
simsk2u4 = sim_pearson_with_k2(mldata,4);
simsk2u5 = sim_pearson_with_k2(mldata,5);


fprintf('User sim\n')
fprintf('User Similarity between 1 and 2:  %.2f\n',simsk2u(2,5))
fprintf('User Similarity between 1 and 4:  %.2f\n',simsk2u(2,7))

% sims2 = corr(mldata', 'rows', 'pairwise');
% fprintf('Item sim\n')
% 
% fprintf('Item Similarity between 1 and 2:  %.2f\n',sims2(7,2))
% fprintf('Item Similarity between 1 and 4:  %.2f\n',sims2(1,4))

%% getting one's sim score

kevin_corrs = simsk2u(4,:);
 kevin_corrs(isnan(kevin_corrs)) = 0;
%  x(x~=x) = 0;
% [ngh_sim, ngh_idx] = sort(kevin_corrs,'descend'); %sorting
[ngh_sim, ngh_idx] = sort(kevin_corrs,'descend');
% ngh_sim = flipud(ngh_sim);
% ngh_idx = flipud(ngh_idx);
% validIndices = ~isnan(ngh_sim);
% ngh_sim = ngh_sim(validIndices);
% kk = zeros(size(kevin_corrs));
% myVec = flipud(myVec)

ngh_sim
ngh_idx

%%
kevin_mu = nanmean(mldata(4,:));           % Kevin's average rating
ngh_corr(4:end) = [];                       % drop non-neighbors
ngh_idx(4:end) = [];                        % drop non-neighbors
ngh_mu = nanmean(mldata(ngh_idx,:),1);       % neighbor average ratings
Predicted = nan(length(ITEM),1);          % initialize an accumulator

for i = 1:length(ITEM)                    % loop over movies
    ngh_r = mldata(ngh_idx,i);             % neighbor ratings for the movie
    isRated = ~isnan(ngh_r);                % only use neighbors who rated
    meanCentered =...                       % mean centered weighted average
        (ngh_r(isRated)' - ngh_mu(isRated)) * ngh_corr(isRated)'...
        / sum(ngh_corr(isRated));
    Predicted(i) = kevin_mu + meanCentered; % add Kevin's average
end

%  Actual = mldata(4,:);                      % Kevin's actual ratings
% table(Actual, Predicted,'RowNames',ITEM)  % compare them to predicted

% fprintf('Predicted rating for "%s": %.d\n',movies{3},round(Predicted(3)))

%% prediction for single user with topN/knn
kevin_mu = nanmean(mldata(4,:));
ngh_sim(4:end) = [];                       % drop non-neighbors
ngh_idx(4:end) = [];                        % drop non-neighbors
ngh_mu = nanmean(mldata(ngh_idx,:),1);
Predicted = zeros(length(ITEM),1);

for i = 1:length(ITEM)                    % loop over movies
    ngh_r = mldata(ngh_idx,i);             % neighbor ratings for the movie
%      isRated = any(ngh_r);                % only use neighbors who rated
       isRated = ~isnan(ngh_r); 
    meanCentered =...                       % mean centered weighted average
        (ngh_r(isRated)' - ngh_mu(isRated)) * ngh_sim(isRated)'...
        / sum(ngh_sim(isRated));
    Predicted(i) =  kevin_mu + meanCentered; % add Kevin's average
end


%%
Actual = mldata(4,:);
RMSE = sqrt(nanmean((Predicted' - Actual).^2))
%%
clear all
clc
