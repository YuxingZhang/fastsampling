
for i = 1:100
    Y = GibbsSampler(X, Y);
    sum(Y == L)
end