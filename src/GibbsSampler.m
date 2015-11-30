function Y = GibbsSampler(X, Y)
% X: D x N
% Y: 1 x N
% N is the total number of data
N = 100;
K = 2;

% X_k is all data in class k
% N_k is a counter for the number of data in each class
% S_k is the sum x_n * x_n^T for each class
% Xsum_k is the sum of data in each class
N_k = [];
S_k = [];
Xsum_k = [];
for k = 1: K
    idx = (Y == k);
    X_k = X(:, idx);
    N_k = [N_k size(X_k, 2)];
    S_k = [S_k X_k * X_k'];
    Xsum_k = [Xsum_k sum(X_k, 2)];
end

%% now build an alias table using the counters above
%% L: 2*c matrix, H: 2*c matrix, A: K bins
L = [];
H = [];
A = [];
l = K;
alpha = K;

%% P is the probability distribution on the K classes
p = (N_k - 1 + (alpha / K)) / (N + alpha - 1);
% normalize p
p = p / sum(p);
for i = 1: l
	pair = [i; p(i)];
	if p(i) <= 1.0 / l
		L = [L pair];
	else 
		H = [H pair];
	end
end

while size(L, 2) > 0
	i = L(1, 1);
    pi = L(2, 1);
	L(:, 1) = []; %% no error checked
	h = H(1, 1);
    ph = H(2, 1);
	H(:, 1) = [];
	bin = [i; h; pi];
	A = [A bin];
	pair = [h; (ph - (1.0 / l - pi))];
	if (ph - (1.0 / l - pi)) > 1.0 / l
		H = [H pair];
	else 
		L = [L pair];
	end
end

%% Now we have the table
k_0 = 5;
v_0 = 5;
m_0 = [0, 0]';
S_0 = [1, 0; 0, 1];
D = 2;

for j = 1: N
    x_n = X(:, j);

    %% Generate a proposal new_k from q(y)
    %new_k = (rand >= 0.5) + 1;
    bin = randi(l);
	i = A(1, bin);
    h = A(2, bin);
    p = A(2, bin);
	if l * p > rand
		new_k = h;
	else
		new_k = i;
	end
    
    % Evaluate the new_k
    new_n = N_k(1, new_k);
    new_k_n = k_0 + new_n;
    new_v_n = v_0 + new_n;
    new_m_n = (k_0 * m_0 + Xsum_k(:, new_k)) / new_k_n;
    new_S_n = S_0 + S_k(:, new_k * D-1:new_k * D) - k_0 * (m_0 * m_0') - new_k_n * (new_m_n * new_m_n');
    
    new_p = GIW(new_k_n, new_v_n, new_S_n, x_n);

    % Evaluate the original k
    k = Y(1, j);
    n = N_k(1, k);
    k_n = k_0 + n;
    v_n = v_0 + n;
    m_n = (k_0 * m_0 + Xsum_k(:, k)) / k_n;
    S_n = S_0 + S_k(:, k * D-1:k * D) - k_0 * (m_0 * m_0') - k_n * (m_n * m_n');
    
    p = GIW(k_n, v_n, S_n, x_n);

    % Evaluate r = min{1,p(yi = new_k | rest) / p(yi = k|rest) * q(k) / q(new_k)}
    r = min([new_p / p, 1]);
    
    % Generate s = Uniform(0; 1) and decide accept or not
    s = rand;
    if s < r
        Y(1, j) = new_k;
        %% TODO updating the counts
        N_k(1, k) = N_k(1, k) - 1;
        N_k(1, new_k) = N_k(1, k) + 1;
    end
end
