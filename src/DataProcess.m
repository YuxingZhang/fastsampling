load '../test.label'
L = test;
load '../test.data'

dim = max(test);
X = zeros(dim(1,1), dim(1,2));
for i = 1 : size(test,1)
    X(test(i, 1), test(i, 2)) = test(i, 3);
end