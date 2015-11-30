load '../train.label'
L = train';
load '../train.data'
data = train;

dim = max(data);
X = zeros(10000, dim(1,1));
for i = 1 : size(data,1)
    if data(i, 2) <= 10000
        X(data(i, 2), data(i, 1)) = data(i, 3);
    end
end