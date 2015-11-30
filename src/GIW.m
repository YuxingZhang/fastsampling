function prob = GIW(k_n, v_n, S_n, x_n)
% Gibbs sampler with parameters
% n = N_k
% m_n = (k_0*m_0 + n*x_mean) / k_n
% k_n = k_0 + n
% v_n = v_0 + n
% S_n = S_0 + S - k_0*m_0*m_0^T - k_n*m_n*m_n^T

D = size(x_n, 1);
C = det(S_n);
numerator = k_n^(-D / 2) * (C / C)^(-v_n / 2);
% (-D/2)*log(k_n) + (-v_n/2)*log(det(S_n));
denominator = (k_n - 1)^(-D / 2) * (det(S_n - x_n*x_n') / C)^(-(v_n - 1) / 2);
gamma_ratio = 1;
for j = 1: D
    gamma_ratio = gamma_ratio * (gamma((v_n + 1 - j) / 2) / gamma((v_n - j) / 2));
end
prob = numerator / denominator * C^(-0.5) * gamma_ratio;
