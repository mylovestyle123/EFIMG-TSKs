function [Result_train,Result_test,W_train,result_x] = prediction_base(W_train, W_test, p, N, N_test, ensemble_size,Temp1, Temp2)

format long

result_x = zeros(N,1,ensemble_size);
result_y = zeros(N_test,1,ensemble_size);
Result_train = zeros(N,ensemble_size);
Result_test = zeros(N_test,ensemble_size);

for i = 1:ensemble_size   
   for j = 1:N
       W_train(j,:,i) = W_train(j,:,i)/(sum(Temp1(j,:,i)));
       result_x(j,1,i) = W_train(j,:,i)*p;
   end
end

for i = 1:ensemble_size
    Result_train(:,i) =  result_x(:,1,i);
end

for i = 1:ensemble_size   
   for j = 1:N_test
       W_test(j,:,i) = W_test(j,:,i)/(sum(Temp2(j,:,i)));
       result_y(j,1,i) = W_test(j,:,i)*p;
   end
end

for i = 1:ensemble_size
    Result_test(:,i) =  result_y(:,1,i);
end

end