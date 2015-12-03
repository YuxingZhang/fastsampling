function Y = Sample(X, Y, K)
N = size(X, 2);
D = size(X, 1);

N_k = [];
S_k = zeros(K, D, D);
Sum_k = [];
for i = 1: K
	X_k = X(:, Y == i);
	N_k = [N_k size(X_k, 2)];
	S_k(i, :, :) = X_k * X_k';
	Sum_k = [Sum_k sum(X_k, 2)];
end

alpha = K;
for i = 1: N
	x_i = X(:, i);

	p = zeros(1, K);
    for k = 1: K
        kappa_0 = 0.01;
        v_0 = D + 2;
        m_0 = Sum_k(:, k) / N_k(1, k);
        X_k = X(:, Y == k);
        S_0 = diag(diag((X_k - repmat(m_0, [1, size(X_k, 2)])) * (X_k - repmat(m_0, [1, size(X_k, 2)]))')) / N_k(1, k);
		if Y(1, i) ~= k
			n = N_k(1, k) + 1;
			kappa_n = kappa_0 + n;
			v_n = v_0 + n;
			m_n = (kappa_0 * m_0 + Sum_k(:, k) + x_i) / kappa_n;
			S_ni = S_0 + squeeze(S_k(k, :, :)) + x_i * x_i' + kappa_0 * (m_0 * m_0') - kappa_n * (m_n * m_n');
            m_n = (kappa_0 * m_0 + Sum_k(:, k)) / (kappa_n - 1);
            S_n = S_0 + squeeze(S_k(k, :, :)) + kappa_0 * (m_0 * m_0') - (kappa_n - 1) * (m_n * m_n');
		else
			n = N_k(1, k);
			kappa_n = kappa_0 + n;
			v_n = v_0 + n;
			m_n = (kappa_0 * m_0 + Sum_k(:, k)) / kappa_n;
			S_ni = S_0 + squeeze(S_k(k, :, :)) + kappa_0 * (m_0 * m_0') - kappa_n * (m_n * m_n');
			m_n = (kappa_0 * m_0 + Sum_k(:, k) - x_i) / (kappa_n - 1);
			S_n = S_0 + squeeze(S_k(k, :, :)) - x_i * x_i' + kappa_0 * (m_0 * m_0') - (kappa_n - 1) * (m_n * m_n');
		end
		p(1, k) = GIW(kappa_n, v_n, S_ni, S_n, x_i) * (n - 1 + (alpha / K)) / (N + alpha - 1);
    end
    p = p / sum(p);

	tmp = rand;
    for k = 1: size(p, 2)
		tmp = tmp - p(1, k);
		if tmp <= 0
			new_k = k;
			break;
		end
    end
    
    k = Y(1, i);
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
