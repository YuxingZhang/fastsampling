function Y = GibbsSampler_o(X, Y, K)
N = size(X, 2);
D = size(X, 1);

% t = cputime
N_k = [];
S_k = zeros(K, D, D);
Sum_k = [];
for i = 1: K
	X_k = X(:, Y == i);
	N_k = [N_k size(X_k, 2)];
	S_k(i, :, :) = X_k * X_k';
	Sum_k = [Sum_k sum(X_k, 2)];
end

% now build an alias table using the counters above
% L: 2*c matrix, H: 2*c matrix, A: K bins
L = [];
H = [];
A = zeros(3, K);
l = K;
alpha = K;

% P is the probability distribution on the K classes
% normalize p
p = (N_k - 1 + (alpha / K)) / (N + alpha - 1);
p = p / sum(p);
for i = 1: l
	pair = [i; p(i)];
	if p(i) <= 1.0 / l
		L = [L pair];
	else 
		H = [H pair];
	end
end

k = 1;
while size(L, 2) > 0
	i = L(1, 1);
    pi = L(2, 1);
	L(:, 1) = []; %% no error checked
    if size(H, 2) > 0
        h = H(1, 1);
        ph = H(2, 1);
        H(:, 1) = [];
    else
        h = i;
        ph = 1;
    end
	bin = [i; h; pi];
	A(:, k) = bin;
	pair = [h; (ph - (1.0 / l - pi))];
	if (ph - (1.0 / l - pi)) > 1.0 / l
		H = [H pair];
	else 
		L = [L pair];
    end
    k = k + 1;
end

while k <= K && size(H, 2) > 0
    i = H(1, 1);
    ph = H(2, 1);
    h = i;
    bin = [i; h; ph];
    A(:, k) = bin;
	H(:, 1) = [];
    k = k + 1;
end

for i = 1: N
    x_i = X(:, i);
    
    % Generate a proposal new_k from q(y)
    bin = randi(l);
    tmp = A(:, bin);
    if l * tmp(3, 1) > rand %% if 
        new_k = tmp(2, 1); %% use h
    else
        new_k = tmp(1, 1); %% use i
    end
    			
    k = Y(1, i);
    
    if new_k == k
                kappa_0 = 0.01;
        v_0 = D + 2;
        m_0 = Sum_k(:, k) / N_k(1, k);
        X_k = X(:, Y == k);
        S_0 = diag(diag((X_k - repmat(m_0, [1, size(X_k, 2)])) * (X_k - repmat(m_0, [1, size(X_k, 2)]))')) / N_k(1, k);
        
        n = N_k(1, k);
		kappa_n = kappa_0 + n;
		v_n = v_0 + n;
		m_n = (kappa_0 * m_0 + Sum_k(:, k)) / kappa_n;
		S_ni = S_0 + squeeze(S_k(k, :, :)) + kappa_0 * (m_0 * m_0') - kappa_n * (m_n * m_n');
		m_n = (kappa_0 * m_0 + Sum_k(:, k) - x_i) / (kappa_n - 1);
		S_n = S_0 + squeeze(S_k(k, :, :)) - x_i * x_i' + kappa_0 * (m_0 * m_0') - (kappa_n - 1) * (m_n * m_n');

        p = GIW(kappa_n, v_n, S_ni, S_n, x_i) * (n - 1 + (alpha / K)) / (N + alpha - 1);
        
        k = new_k;
        kappa_0 = 0.01;
        v_0 = D + 2;
        m_0 = Sum_k(:, k) / N_k(1, k);
        X_k = X(:, Y == k);
        S_0 = diag(diag((X_k - repmat(m_0, [1, size(X_k, 2)])) * (X_k - repmat(m_0, [1, size(X_k, 2)]))')) / N_k(1, k);
        
        n = N_k(1, k);
		kappa_n = kappa_0 + n;
		v_n = v_0 + n;
		m_n = (kappa_0 * m_0 + Sum_k(:, k)) / kappa_n;
		S_ni = S_0 + squeeze(S_k(k, :, :)) + kappa_0 * (m_0 * m_0') - kappa_n * (m_n * m_n');
		m_n = (kappa_0 * m_0 + Sum_k(:, k) - x_i) / (kappa_n - 1);
		S_n = S_0 + squeeze(S_k(k, :, :)) - x_i * x_i' + kappa_0 * (m_0 * m_0') - (kappa_n - 1) * (m_n * m_n');

        new_p = GIW(kappa_n, v_n, S_ni, S_n, x_i) * (n - 1 + (alpha / K)) / (N + alpha - 1);
        k = Y(1, i);
        
        r = 1;
    else
        kappa_0 = 0.01;
        v_0 = D + 2;
        m_0 = Sum_k(:, k) / N_k(1, k);
        X_k = X(:, Y == k);
        S_0 = diag(diag((X_k - repmat(m_0, [1, size(X_k, 2)])) * (X_k - repmat(m_0, [1, size(X_k, 2)]))')) / N_k(1, k);
        
        n = N_k(1, k);
		kappa_n = kappa_0 + n;
		v_n = v_0 + n;
		m_n = (kappa_0 * m_0 + Sum_k(:, k)) / kappa_n;
		S_ni = S_0 + squeeze(S_k(k, :, :)) + kappa_0 * (m_0 * m_0') - kappa_n * (m_n * m_n');
		m_n = (kappa_0 * m_0 + Sum_k(:, k) - x_i) / (kappa_n - 1);
		S_n = S_0 + squeeze(S_k(k, :, :)) - x_i * x_i' + kappa_0 * (m_0 * m_0') - (kappa_n - 1) * (m_n * m_n');

        p = GIW(kappa_n, v_n, S_ni, S_n, x_i);
        
        k = new_k;
        kappa_0 = 0.01;
        v_0 = D + 2;
        m_0 = Sum_k(:, k) / N_k(1, k);
        X_k = X(:, Y == k);
        S_0 = diag(diag((X_k - repmat(m_0, [1, size(X_k, 2)])) * (X_k - repmat(m_0, [1, size(X_k, 2)]))')) / N_k(1, k);

        n = N_k(1, k) + 1;
		kappa_n = kappa_0 + n;
		v_n = v_0 + n;
		m_n = (kappa_0 * m_0 + Sum_k(:, k) + x_i) / kappa_n;
		S_ni = S_0 + squeeze(S_k(k, :, :)) + x_i * x_i' + kappa_0 * (m_0 * m_0') - kappa_n * (m_n * m_n');
        m_n = (kappa_0 * m_0 + Sum_k(:, k)) / (kappa_n - 1);
        S_n = S_0 + squeeze(S_k(k, :, :)) + kappa_0 * (m_0 * m_0') - (kappa_n - 1) * (m_n * m_n');
		
        new_p = GIW(kappa_n, v_n, S_ni, S_n, x_i);
        k = Y(1, i);
        
		% Evaluate r = min{1,p(yi = new_k | rest) / p(yi = k|rest) * q(k) / q(new_k)}
        r = min([new_p / p, 1]);        
    end
    % Generate s = Uniform(0; 1) and decide accept or not
    s = rand;
    
    if s < r && new_k ~= k
        Y(1, i) = new_k;
        X_k = zeros(1, D, D);
        X_k(1, :, :) = x_i * x_i';
        S_k(k , :, :) = S_k(k, :, :) - X_k;
        S_k(new_k, :, :) = S_k(new_k, :, :) + X_k;
        Sum_k(:, k) = Sum_k(:, k) - x_i;
        Sum_k(:, new_k) = Sum_k(:, new_k) + x_i;
        N_k(1, k) = N_k(1, k) - 1;
        N_k(1, new_k) = N_k(1, new_k) + 1;
    end
end

% cputime - t
