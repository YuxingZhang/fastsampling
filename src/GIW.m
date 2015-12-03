function prob = GIW(k_n, v_n, S_ni, S_n, x_n)
% Gibbs sampler with parameters
% n = N_k
% m_n = (k_0*m_0 + n*x_mean) / k_n
% k_n = k_0 + n
% v_n = v_0 + n
% S_n = S_0 + S - k_0*m_0*m_0^T - k_n*m_n*m_n^T

D = size(x_n, 1);
C = det(S_ni);
numerator = k_n^(-D / 2) * (C / C)^(-v_n / 2);
denominator = (k_n - 1)^(-D / 2) * (det(S_n) / C)^(-(v_n - 1) / 2);
gamma_ratio = 0;
for j = 1: D
    gamma_ratio = gamma_ratio + gammaln((v_n + 1 - j) / 2) - gammaln((v_n - j) / 2);
end
prob = numerator / denominator * C^(-0.5) * exp(gamma_ratio);

% numerator = (-D / 2) * log(kappa_n) + (-v_n / 2) * log(C / C);
% denominator = (-D / 2) * log(kappa_n - 1) + (-(v_n - 1) / 2) * log(abs(det(S_n - x_n*x_n') / C));
% gamma_ratio = 0;
% for j = 1: D
%     gamma_ratio = gamma_ratio + gammaln((v_n + 1 - j) / 2) - gammaln((v_n - j) / 2);
% end
% prob = numerator - denominator + log(C^(-0.5)) + gamma_ratio;
