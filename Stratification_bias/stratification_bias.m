m = 200;
n = 100;

h2 = 0.3;
sigm_s = 0.3;
rho = 0.3;

% Create the covariance matrix Sigm
Sigm = rho.^abs(bsxfun(@minus, (1:m)', 1:m));

% Perform eigen decomposition
[eigenvectors, eigenvalues] = eig(Sigm);

Sigm_half = eigenvectors * diag(diag(eigenvalues).^(1/2)) * eigenvectors';

n_sim = 10^3;
res_reg_free = NaN(n_sim, 1);
do_scale = false;

F_ST_ = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];
n_F_ST = length(F_ST_);

res_m_sim = NaN(n_F_ST, 1);
res_s_sim = NaN(n_F_ST, 1);
res_theory = NaN(n_F_ST, 1);

for v = 1:n_F_ST
    disp(v);
    F_ST = F_ST_(v);
    f_vec = randn(m, 1) * sqrt(F_ST);

    for r = 1:n_sim
        X = NaN(n, m);
        Z = reshape(randn(n*m, 1), [n, m]);

        for i = 1:n
            if i <= n/2
                X(i, :) = ((Sigm_half * Z(i, :)')' - f_vec') / sqrt(1 + f_vec'.^2);
            else
                X(i, :) = ((Sigm_half * Z(i, :)')' + f_vec') / sqrt(1 + f_vec'.^2);
            end
        end
        bet = randn(m, 1) * sqrt(h2/m);
        xi = [-sigm_s * ones(n/2, 1); sigm_s * ones(n/2, 1)];
        Y = X * bet + xi + randn(n, 1) * sqrt(1 - h2 - sigm_s^2);

        V_y = var(Y);

        if do_scale
            X = X - mean(X);
            Y = Y / sqrt(V_y);

            u2 = (X' * Y).^2 / (n - 1);
            R = (X' * X) / (n - 1);
        else
            u2 = (X' * Y).^2 / n;
            R = (X' * X) / n;
        end

        hatell = sum(R.^2, 2);

        res_reg_free(r) = mean(u2 .* (n/m*hatell - mean(n/m*hatell))) / mean((n/m*hatell - mean(n/m*hatell)).^2);
    end

    res_m_sim(v) = mean(res_reg_free) - h2;
    res_s_sim(v) = std(res_reg_free) / sqrt(n_sim);
    C_f = mean(f_vec.^2 ./ (1 + f_vec.^2));
    res_theory(v) = sigm_s^2 / C_f;
end

res_m_sim
res_theory

% Plot the results
figure;
plot(F_ST_, res_m_sim, 'LineWidth', 2);
xlabel('F_{ST}');
ylabel('bias');
title('m=200,n=100,Stand.');
ylim([min(res_m_sim - 2*res_s_sim), max(res_m_sim + 2*res_s_sim)]);

hold on;
for j = 1:n_F_ST
    line([F_ST_(j), F_ST_(j)], [res_m_sim(j) - 2*res_s_sim(j), res_m_sim(j) + 2*res_s_sim(j)], 'LineStyle', '--');
end

plot(F_ST_, res_theory, 'r');
hold off;
