function [traindata_x, traindata_y,testdata_x,testdata_y] = preprocess(traindata,testdata,N,D)

% random sampling

trainingsubset = zeros(N,D+1);
temp1 = round(rand(N,1)*(N-1))+1; 

for i = 1:N
     trainingsubset(i,:) = traindata(temp1(i,1),:);
end

traindata_x = trainingsubset(:,1:D);
traindata_y = trainingsubset(:,D+1);
testdata_x = testdata(:,1:D);
testdata_y = testdata(:,D+1);

end