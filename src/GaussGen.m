close all;

w1 = rand;
mu1 = [0, 10];
mu2 = [10, 0];
Sigma1 = [1, 0; 0, 1];
Sigma2 = [2, 0; 0, 6];
X = [];
L = [];
hold on;
for i = 1:100
    w = rand;
    if w >= w1
        X = [X mvnrnd(mu1, Sigma1)'];
        L = [L 1];
    else
        X = [X mvnrnd(mu2, Sigma2)'];
        L = [L 2];
    end
end
scatter(X(1,:), X(2,:));
D = size(X, 1);
N = size(X, 2);
K = 2;
Y = randi([1, K], 1, N);
    
