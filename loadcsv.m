%% load csv into matrix
% filename = '/Users/xahiru/Documents/MATLAB/ml-latest-small/ratings.csv';
filename = '/Users/xahiru/Documents/MATLAB/ml-latest-small/ratings2.csv';
DATA = csvread(filename,1,0);
save('DATA');

clear filename;

%% sort data
% C is the user
% B is items
[USER,ia,ic] = unique(DATA(:,1),'rows');
[ITEM,ik,ib] = unique(DATA(:,2),'rows');

%ia is the index of unique id users, ic is the index of all ids the users
%ik is the index of unique id moives, ic is the index of all ids the movies

% create the zero array to hold data
mldata=zeros(size(ia,1),size(ik,1));

% go through all the rows to create the data matrix
for j=1:size(DATA(:,1))
    user= find(USER(:,1)== DATA(j,1));
    movie= find(ITEM(:,1)== DATA(j,2));
    rating= DATA(j,3);
    mldata(user,movie)=rating;
    
end

clear rating movie user i
%% Using the MovieLens Dataset - build the item-similarity dataset
%
% This will take a while (10 minutes or so) so go take a break. 

 disp('Building item-similarity dataset...');
 tStart = tic;
% 
% itemsim=calculateSimilarItems(mldata,50,@sim_distance);
  itemsim=calculateSimilarItems(mldata,50,@sim_pearson);
%   itemsim_p_k = calculateSimilarItems(mldata,50,@sim_pearson_with_k);
 
 tEnd = toc(tStart);
 fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
 
%  save('itemsim_p_k');
 save('itemsim');
 
 clear tEnd tStart
 %% Recommending Items (Page 17)
%
% Similarity scores can be used to make recommendations. You can find who
% has similar tastes to given person and find out what other movies those
% people have seen that the given person has not seen. But we need to use
% weighted average to make the result consistent. "recommendations for
% Toby.xls" shows how this algorithm works in Excel.

% get user-based recommendations for Toby
  recom = getRecommendations(mldata,7,@sim_pearson_with_k,2);
%  
%  temp = getRecommendations(mldata,7,@sim_pearson,2);


disp('Recommend movies for 7 - user-based filtering with sim_pearson_with_k')
disp(recom)

clear i recom ;

%% Recommendation for all

% temp = getRecommendations(mldata,7,@sim_pearson_with_k,2);



tempCell = cell((USER:1),3);


 for j=1:size(USER(:,1))
  
     recs =  getRecommendations(mldata,j,@sim_pearson_with_k,2);
     
     
  if ~isempty(recs)
     tempCell{j,1}=recs(:,1);
     tempCell{j,2}=recs(:,2);
     tempCell{j,3} = USER(j,:);
     
     Predicted(USER(j,:),recs(:,2)) = recs(:,1);
  else
      tempCell{j,1} = NaN;
  end   
 

  
  
 end
 
 clear  i
 
 %% Using the MovieLens Dataset - get item-based recommendations

% get item-based recommendations for User ID = 87
% temp=getRecommendedItems(prefs,itemsim,87);
temp=getRecommendedItems(mldata,itemsim,7);


tempCell = cell(size(temp,1),2);

for j=1:size(temp,1)
    tempCell{j,1}=temp(j,1);
    tempCell{j,2}=temp(j,2);
%     tempCell{i,2}=ITEM(1,temp(i,2));
end

disp('Recommend movies for User [ID = 87] - item-based filtering')
disp(tempCell(1:10,:))


clear i 

%% predictALl items for users

Predicteditems = zeros(size(copymldata));
for i=1:size(USER(:,1))
temp=getRecommendedItems(copymldata,itemsim_p_k,i);


   if ~isempty(temp)
     
      Predicteditems(USER(i,:),temp(:,2)) = temp(:,1);
   else
       i
   end      
        
  temp
        

end
clear i j tempCell temp

%% mean centering

[m,n] = size(mldata);

meanRatings = sum(mldata,2)./sum(mldata~=0,2);
[i,j,v] = find(mldata);
meanCentered = sparse(i,j,v - meanRatings(i),m,n);


itemsim_p_k_mc =calculateSimilarItems(meanCentered,50,@sim_pearson_with_k);
save('itemsim_p_k_mc');

clear i m n j v meanRatings
%% NNMatrix FActoriztion
 k = 4;
 [W,H] = nnmf(mldata,k);
 
 
 %% 10 Fold crossval
load('DATA');
CVO = cvpartition(mldata(1,:),'k',10);

err = zeros(CVO.NumTestSets,1);
for j = 1:CVO.NumTestSets
    trIdx = CVO.training(j);
    teIdx = CVO.test(j);
    ytest = classify(meas(teIdx,:),meas(trIdx,:),...
		 ratings(trIdx,:));
    err(j) = sum(~strcmp(ytest,mldata(teIdx)));
end
cvErr = sum(err)/sum(CVO.TestSize);

 %% Topmatch test
%  scores = topMatches(mldata',1,5,@sim_pearson_with_k,2)

% scores = sim_pearson(mldata,x1,x2,0)
itemsim =calculateSimilarItems(mldata,50,@sim_pearson,1)

x1 = 1;
x2 = 2;
    c1=find(mldata(x1,:)>0);
    c2=find(mldata(x2,:)>0);
    cc=intersect(c1,c2);
    shared_items(1,:)=mldata(x1,cc);
    shared_items(2,:)=mldata(x2,cc);
    
    % get the size of shared_items
    n=size(shared_items,2);
    n

 clear c1 c2 cc
 
 %% test area
mat2 = zeros(size(mldata));

copymldata = mldata;

copymldata(2, 3:20) = 0;

%  mat2(:, 20:40) = mldata(:, 20:40)

mat2 = ;
 

 %% clear
clear all
clear
clc
