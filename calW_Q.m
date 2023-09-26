function [w, q,temp3] = calW_Q(mu, sigma, variance, X)
format long

% mu: (Nc,D) 
% sigama: (D,D,Nc)
% X: (N,D)
N = size(X,1);
d = size(X,2);
Nc = size(mu, 1);
D = size(mu, 2);

w = zeros(N, Nc);
q = zeros(N, Nc);
temp3 = zeros(N,Nc);
% temp4 = zeros(N,1);

% calculate every coefficient of w
for i = 1:N 
    x = X(i, :);  % x: (1,d)
    for j = 1:Nc
        mu_cx = mu(j, 1:d);  % mu_cx: (1,d)
        sigma_cxy = sigma(1:d, D, j);  % sigma_cxy: (d,1)
        sigma_cyy = sigma(D, D, j);  % sigma_cyy: (1,1)
        
        % size(x)
        % size(mu_cx)
        a = x - mu_cx;  % a: (1,d)
        b = sigma_cxy/sigma_cyy;  % b: (d,1)
        
        temp1 = (mu(j, D) - dot(a, b));
        % temp1
        
        temp2 = 1;
        % calculate temp2 according to Eq.(5)
        for k = 1:d
            xd = x(1, k);
            temp2 = normpdf(xd, mu(j, k), sqrt(variance(j, k))) * temp2;  
        end
        
        % calculate W and Q according to Eqs.(8) and (9)
        w(i, j) = temp1 * temp2;  
        q(i, j) = (temp1 * temp1 + 1/sigma_cyy) * temp2;

        temp3(i,j) = temp2/Nc;  
    end
    

end
