for i = 1:10
    Y = GibbsSampler(X, Y);
    sum(Y == L)
end