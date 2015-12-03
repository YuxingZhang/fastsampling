for i = 1 : 100
    i
%     Y = Sample(X, Y);
%     Y = GibbsSampler(X, Y);
     Y = GibbsSampler(X, Y, 10);
end
K = 4;
for i = 1 : K
    idx = (Y == i);
    X1 = X(:, idx);
    scatter(X1(1, :), X1(2,:));
end