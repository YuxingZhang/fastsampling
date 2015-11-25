%% The main program for the fast sampling Gaussian mixture model


%% Input format: D×N matrix, where D is the dimension of the data, 
%% N is the number of data points we have
filename = 'InputDataSet.dat';
X = csvread(filename);
D = size(X, 1);
N = size(X, 2);

%% K is the number of classes (number of mixture components) which 
%% is specified by the user.
K = 10;

%% random sample the class labels initially
Y = randi([1, K], 1, N);

%% T is the number of iterations we sample
T = 100;

for t = 1 : T
    
end









