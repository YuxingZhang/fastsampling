function prob = GIW(k_n, v_n, S_n, x_n)
% Gibbs sampler with parameters
% n = N_k
% m_n = (k_0*m_0 + n*x_mean) / k_n
% k_n = k_0 + n
% v_n = v_0 + n
% S_n = S_0 + S - k_0*m_0*m_0^T - k_n*m_n*m_n^T

D = size(x_n, 1);
numerator = k_n^(-D/2) * det(S_n)^(-v_n/2);
denominator = (k_n-1)^(-D/2) * det(S_n - x_n*x_n')^(-(v_n-1)/2);
for j = 1:D
    numerator = numerator * gamma((v_n+1-j)/2);
    denominator = denominator * gamma((v_n-j)/2);
end
prob = numerator / denominator;