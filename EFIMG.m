%% Load data
clear;
close all;
% training data
[traindata] = xlsread('ccpp_train.xlsx',1,'A1£ºE6698');
train_x = traindata(:,1:4);
train_y = traindata(:,5);
N = size(train_x,1); %training samples
D = size(train_x,2);
% testing data
[testdata] = xlsread('ccpp_test.xlsx',1,'A1£ºE2870');
testy = testdata(:,5);
N_test = size(testdata,1); %testing samples
D_test = size(testdata,2)-1;

%% Preprocess for base FIMG-TSKs
% % normalizeData--the data we used have already been normalized.
% traindata = [normalizeData(train_x) train_y];
% testdata = [normalizeData(testdata(:,1:4)) testy];

% hyperparameters--ensemble size
ensemble_size = 6;

ensemble_train_record = zeros(N, D+1, ensemble_size);
ensemble_test_record = zeros(N_test, D_test+1, ensemble_size);

% function "preprocess"
for i = 1 : ensemble_size
    [traindata_x, traindata_y,testdata_x,testdata_y] =  preprocess(traindata,testdata,N,D);
    ensemble_train_record(1:N,:,i) = [traindata_x, traindata_y];
    ensemble_test_record(1:N_test,:,i) = [testdata_x, testdata_y];
end

%% The related parameters for all base FIMG-TSKs by FCM
% hyperparameters--the clusters (number of rules) for each base learner
% Here we uniformly set the number of rules for each learner as Nc for simplicity.
Nc = zeros(ensemble_size,1);
for i = 1 : ensemble_size
    Nc(i,1) = 5;
end

mu_record = zeros(Nc(1,1), D+1, ensemble_size);
sigma_record = zeros(D+1, D+1, Nc(1,1), ensemble_size);
variance_record = zeros(Nc(1,1), D+1, ensemble_size);

% function "myFCM"
for i = 1:ensemble_size
   
    [mu, sigma, variance, covariance, Xrecord, centers, U, row, num] = myFCM(ensemble_train_record(1:N,1:D,i), ensemble_train_record(1:N,D+1,i), Nc(1,1));
    mu_record(1:Nc(1,1),:,i) = mu;
    sigma_record(:,:,:,i) = sigma;
    variance_record(:,:,i) = variance;
    
end

%% Calculate W and Q according to Eqs.(8) and (9)

W_train = zeros(N,Nc(1,1),ensemble_size);
Q_train = zeros(N,Nc(1,1),ensemble_size);
W_test = zeros(N_test,Nc(1,1),ensemble_size);
Q_test = zeros(N_test,Nc(1,1),ensemble_size);
Temp1 = zeros(N,Nc(1,1),ensemble_size);
Temp2 = zeros(N_test,Nc(1,1),ensemble_size);

% function "calW_Q"
for i = 1:ensemble_size
    [w_train, q_train,temp1] = calW_Q(mu_record(1:Nc(1,1),:,i), sigma_record(:,:,:,i), variance_record(:,:,i), ensemble_train_record(1:N,1:D,i));
    [w_test, q_test,temp2] = calW_Q(mu_record(1:Nc(1,1),:,i), sigma_record(:,:,:,i), variance_record(:,:,i), ensemble_test_record(1:N_test,1:D,i));
    W_train(:,:,i) = w_train;
    Q_train(:,:,i) = q_train;
    W_test(:,:,i) = w_test;
    Q_test(:,:,i) = q_test;
    Temp1(:,:,i) = temp1;
    Temp2(:,:,i) = temp2;
end

%% Calculate P

% The functions "cal_Z","cal_D","calW_W" and "cal_p" are saved in our personal computers.
% If you need, please contact us.

% p = zeros(sum(Nc,1),1);

% %%% define matrix H according to Eq.(18)
% H = zeros(sum(Nc,1),sum(Nc,1),ensemble_size);
% k = 0;
% for i = 1: ensemble_size
%     for j = 1:Nc(i,1)
%         H(j+k,j+k,i) = 1; 
%     end
%     k = k+Nc(i,1);
% end

% %%% define matrix Z according to Eq.(15)
% function "cal_Z"
% k = 1;
% Z = zeros(sum(Nc,1),sum(Nc,1),ensemble_size);
% for i = 1:(ensemble_size-1)
%     for j = (i+1):ensemble_size
%         Z(:,:,k) = cal_Z(mu_record(1:Nc(1,1),:,i), sigma_record(:,:,:,i), variance_record(:,:,i), ensemble_train_record(1:N,1:D,i),mu_record(1:Nc(1,1),:,j), sigma_record(:,:,:,j), variance_record(:,:,j), ensemble_train_record(1:N,1:D,j));
%         k = k+1; 
%     end
% end

% %%% define matrix D according to Eq.(20)
% function "cal_D"
% k = 1;
% D = zeros(sum(Nc,1),sum(Nc,1),ensemble_size);
% for i = 1:(ensemble_size-1)
%     for j = (i+1):ensemble_size
%         D(:,:,k) = cal_D(mu_record(1:Nc(1,1),:,i), sigma_record(:,:,:,i), ensemble_train_record(1:N,1:D,i),mu_record(1:Nc(1,1),:,j), sigma_record(:,:,:,j), ensemble_train_record(1:N,1:D,j));
%         k = k+1; 
%     end
% end

% %%% define W_W_train according to Eq.(16) and Eq.(17)
% function "calW_W"
% W_W_train = zeros(sum(Nc(:,1)),1,ensemble_size-1);
% k = 1;
% for i = 1:(ensemble_size-1)
%     for j = (i+1):ensemble_size
%     [w_w_train] = calW_W(mu_record(1:Nc(1,1),:,i), sigma_record(:,:,:,i), variance_record(:,:,i), ensemble_train_record(1:N,1:D,i),mu_record(1:Nc(1,1),:,j), sigma_record(:,:,:,j), variance_record(:,:,j), ensemble_train_record(1:N,1:D,j));
%     W_W_train(:,:,k) = w_w_train;
%     k = k+1;
%     end
% end

% %%% calculate p
% function "cal_p"
% % given hyperparameters
% alpha = 0.1; betta = 100; gamma = 1e+03; yita = 1.2;

% % calculate p according to Eq.(27)
% [p] = cal_p(W_train, Q_train, Z, D, W_W_train, test_y, alpha, betta, gamma, yita);

%% We set the weight of each rule in each base FIMG-TSK as 1/Nc for simplicity.
p=(ones(Nc(1,1),1))*(1/Nc(1,1));

%% Prediction of EFIMG-TSKs
% prediction of training and testing results of each base FIMG-TSK

% function "prediction_base"
[Result_train,Result_test] = prediction_base(W_train, W_test, p, N, N_test, ensemble_size, Temp1, Temp2);

% testing RMSE of the model on CCPP
result = zeros(N_test,1);
for i = 1:N_test
    result(i,1) = sum(Result_test(i,:))/ensemble_size;
end

RMSE_test = sqrt(((testy-result)'*(testy-result))/N_test)






