function [X, Y, L] = newGaussGen(K, N)
close all;

w = zeros(1, K);
mu = zeros(2, K);
Sigma = zeros(2, 2 * K);
for i = 1:K
    w(1, i) = rand;
    mu(:, i) = unidrnd(500, [1,2]);
    Sigma(:, (2*i-1):2*i) = diag(unidrnd(100, [2,1]));
end
w = w / sum(w);

X = [];
L = [];
hold on;
for i = 1:N
    tmp = rand;
    for k = 1:K
        tmp = tmp - w(1, k);
        if tmp <= 0
            X = [X mvnrnd(mu(:, k), Sigma(:, (2*k-1):2*k))'];
            L = [L k];
            break;
        end
    end
end

scatter(X(1,:), X(2,:));
Y = randi([1, K], 1, N);
    
