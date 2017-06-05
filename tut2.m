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
  mldata=NaN(size(ia,1),size(ik,1));
%  mldata=zeros(size(ia,1),size(ik,1));
% go through all the rows to create the data matrix
for j=1:size(DATA(:,1))
    user= find(USER(:,1)== DATA(j,1));
    movie= find(ITEM(:,1)== DATA(j,2));
    rating= DATA(j,3);
    mldata(user,movie)=rating;
    
end

clear rating movie user i

       

%%  correlation/similarity with k user / items
% sims = corr(mldata, 'rows', 'pairwise');


   simsk2u0 = sim_pearson_with_k2(mldata,0);
   simsk2u1 = sim_pearson_with_k2(mldata,1);
   simsk2u2 = sim_pearson_with_k2(mldata,2);
   simsk2u3 = sim_pearson_with_k2(mldata,3);
   simsk2u4 = sim_pearson_with_k2(mldata,4);
   simsk2u5 = sim_pearson_with_k2(mldata,5);
   simown = sim_pearson_with_k2(mldata,6);

  simsk2u = simsk2u4;
%     q = simsk2u4;

%  simsk2u(logical(eye(size(simsk2u)))) = [];

% simsk2u = reshape(simsk2u,length(q)-1, length(q));

% fprintf('User sim\n')
% fprintf('User Similarity between 1 and 2:  %.2f\n',simsk2u(2,5))
% fprintf('User Similarity between 1 and 4:  %.2f\n',simsk2u(2,7))

% sims2 = corr(mldata', 'rows', 'pairwise');
% fprintf('Item sim\n')
% 
% fprintf('Item Similarity between 1 and 2:  %.2f\n',sims2(7,2))
% fprintf('Item Similarity between 1 and 4:  %.2f\n',sims2(1,4))

% simsk2u = simsk2u - eye(length(USER)); % set self-correlations to 0

%% getting one's sim score

usr = 3;
kevin_corrs = simsk2u(usr,:);
%  kevin_corrs(isnan(kevin_corrs)) = 0;
kevin_corrs(usr) = [];

[ngh_sim, ngh_idx] = sort(kevin_corrs,'descend');


%  x(x~=x) = 0;
% [ngh_sim, ngh_idx] = sort(kevin_corrs,'descend'); %sorting

% ngh_sim = flipud(ngh_sim);
% ngh_idx = flipud(ngh_idx);
% validIndices = ~isnan(ngh_sim);
% ngh_sim = ngh_sim(validIndices);
% kk = zeros(size(kevin_corrs));
% myVec = flipud(myVec)

   ngh_sim
   ngh_idx


%% prediction for single user with topN/knn


kevin_mu = nanmean(mldata(4,:));
ngh_sim(2:end) = [];                       % drop non-neighbors
ngh_idx(2:end) = [];                        % drop non-neighbors

% result   = x(sort(index(1:n)))

ngh_mu = nanmean(mldata(ngh_idx,:),1);
Predicted = nan(length(ITEM),1);

for i = 1:length(ITEM)                    % loop over movies
    ngh_r = mldata(ngh_idx,i);             % neighbor ratings for the movie
%      isRated = any(ngh_r);                % only use neighbors who rated
       isRated = ~isnan(ngh_r); 
    meanCentered =...                       % mean centered weighted average
        (ngh_r(isRated)' - ngh_mu(isRated)) * ngh_sim(isRated)'...
        / sum(ngh_sim(isRated));

    Predicted(i) =  kevin_mu + meanCentered; % add Kevin's average
end


%% testing error
usr = 4;
kevin_corrs = simsk2u(usr,:);
%  kevin_corrs(isnan(kevin_corrs)) = 0;
kevin_corrs(usr) = [];
kevin_mu = nanmean(mldata(usr,:));

[a, b] = sort(kevin_corrs,'descend');

KN = 1;
cello = 1;
ngh_sim = a(1:KN);
ngh_idx = b(1:KN);
ngh_mu = nanmean(mldata(ngh_idx,:));

ngh_r = mldata(ngh_idx,cello);

isRated = ~isnan(ngh_r); 
     meanCentered =...                       % mean centered weighted average
         (ngh_r(isRated)' - ngh_mu(isRated)) * ngh_sim(isRated)'...
         / sum(ngh_sim(isRated));
 

%     if isempty(meanCentered)
%     meanCentered = NaN;
%     end


 Predicted(i) =  kevin_mu + meanCentered;

 %% ii

       BigPridcted = nan(size(mldata));
      KN = 6; %should be less than 6
      L = length(ITEM);

    
      for i=1:length(USER)
 %           ///predictWithKN(i,mldata,KN,simsk2u,L)';
 %          //if (~isnan(QQ))
        BigPridcted(i,:) = predictWithKN(i,mldata,KN,simsk2u,L)';
 %         // end
      end
%%
Actual = mldata(4,:);
RMSE = sqrt(nanmean((Predicted' - Actual).^2))

%%

W = ((BigPridcted-mldata).^2);

bigmean = nanmean(W(:)');

RMSE =  sqrt(bigmean);

%% Saving RMSE

   BigPridcted = nan(size(mldata));
%       KN = 3; %should be less than 6
      L = length(ITEM);

   for j=1:length(USER)-1

        
      for i=1:length(USER)
          
        BigPridcted(i,:) = predictWithKN(i,mldata,j,simsk2u,L)';

        
       end

        W = ((BigPridcted-mldata).^2);

        bigmean = nanmean(W(:)');
        
        RMSE(j,:) =  sqrt(bigmean);

    end
%% plotting results
