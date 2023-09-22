function [ mu, sigma, variance, covariance, Xrecord, centers, U, row, num ] = myFCM(xTrain, yTrain, Nc)

format long

if (size(xTrain,1)~=size(yTrain,1))
     fprintf(2, strcat('\n****************************************************',...
        '\nERROR: Number of features of training and test data are not same!',...
        '\n****************************************************\n'));
    return;
end

N = size(xTrain,1);

% concatenate xTrain and yTrain
X = [xTrain, yTrain];  % X: (N,D+1)
D = size(xTrain,2) + 1;  % D + 1

% apply FCM on X
options = [2,100,1e-5,0]; %default settings in FCM
[centers, U] = fcm(X, Nc,options); 

% centers: (Nc,D) -cluster center for every cluster on every dimension.
% U: (Nc, N)      -membership degree for every sample on every cluster.

%mu = centers; % (Nc,D) for every cluster, there are D dimensions. 
mu = zeros(Nc, D);
sigma = zeros(D, D, Nc);
variance = zeros(Nc, D);
covariance = zeros(D, D, Nc);
Xrecord = zeros(N, D, Nc);

% find the corresponding clusters for N samples
value = max(U);
[row, ~] = find(U == value); % row: (N,1)  col: (N,1)

num = zeros(Nc, 1);  % numbers of samples in Nc clusters.
for i = 1:N
    cluster = row(i,:);
    num(cluster,:) = num(cluster,:) + 1;
end

for i = 1:Nc
    Ni = num(i,1);  % number of samples in this cluster
    tempX = zeros(Ni,D);  % (Ni,D)
    ittr = 1;
    for j = 1:N  % check all samples
        if row(j,:) == i  % if it belongs to ith cluster
            tempX(ittr,:) = X(j,:);
            ittr = ittr + 1;
        end
    end
    % size(tempX)
    Xrecord(1:Ni, :, i) = tempX;

    % get sigma for every tempX
    covtemp = cov(tempX) + 1e-7;
    covariance(:, :, i) = covtemp;

    %size(covariance)
    
    variance(i, :) = diag(covtemp);
    sigma(:, :, i) = inv(covtemp);  % sigma: (D,D)
    
    mu(i, :) = mean(tempX);
end






















