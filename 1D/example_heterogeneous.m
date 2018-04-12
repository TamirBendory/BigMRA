clear;
close all;
clc;

%% Pick parameters and generate signals

% Pick K signals of length L and the size W of the separation window
K = 3;
L = 21;
W = 2*L-1;
X = randn(L, K);  % randn instead of rand here -- to be discussed
X(:, 1) = [ones(ceil(L/2), 1) ; -ones(floor(L/2), 1)];
X(:, 2) = [linspace(0, 2, ceil(L/2))' ; linspace(2, 0, floor(L/2))'];
X(:, 3) = randn(L, 1);

% Pick a noise level
sigma = .0;

% Desired number of occurrences of each signal X(:, k)
m_want = [3e6 2e6 1e6] / 1e3;

% Length of micrograph
n = sum(m_want)*W*10;

fprintf('Micrograph length: %g\n\n\n', n);

%% Pick which correlation coefficients to sample

[list2, list3] = moment_selection(L, 'include biased');

%% Generate the micrograph

T = tic();
[y_clean, m_actual] = generate_clean_micrograph_1D_heterogeneous(X, W, n, m_want);
y_obs = y_clean + sigma*randn(n, 1);
time_to_generate_micrograph = toc(T);
fprintf('Time to generate micrograph: %.2g [s]\n', time_to_generate_micrograph);
SNR = norm(y_clean, 'fro')/norm(y_obs-y_clean, 'fro');
fprintf('   SNR: %.2g\n', SNR);
fprintf('   m_actual/m_want: ');
fprintf(' %.2g', m_actual./m_want);
fprintf('\n');

%% Collect the moments
T = tic();
% [M1, M2, M3] = moments_from_data_no_debias_1D(y_obs, list2, list3);
batch_size = 1e8;
[M1, M2, M3] = moments_from_data_no_debias_1D_batch(y_obs, list2, list3, batch_size);
time_to_compute_moments = toc(T);
fprintf('   Moment computation: %.4g [s]\n', time_to_compute_moments);

moments.M1 = M1 / n;  %%%% !!!!! We normalize by n here
moments.M2 = M2 / n;
moments.M3 = M3 / n;
moments.list2 = list2;
moments.list3 = list3;

clear y_clean y_obs;
save(sprintf('data_example_heterogeneous_%s_n_%d.mat', datestr(now(), 30), n));

%% Optimization

% Here we can choose the length of sought signals in optimization.
% Should be no less than L; typically set to W.
L_optim = 2*L-1;

% X0 = randn(L_optim, K);
% gamma0 = m_actual*L_optim/n; % give true gamma for now
sigma_est = sigma; % irrelevant if biased terms are excluded and if the weights internally do not depend on sigma

[X_est, gamma_est, problem1, stats1] = least_squares_1D_heterogeneous(moments, L_optim, K, sigma_est); %, X0, gamma0(:)); % check if the code fixes gamma to gamma0 or not


%% Display

fprintf('Estimated densities:\n');
disp(gamma_est');
fprintf('True densities:\n');
disp(m_actual*L_optim/n);

X_extended = [X ; zeros(L_optim-L, K)];

figure(1);
clf;

for k1 = 1 : K
    for k2 = 1 : K
        subplot(K, K, (k1-1)*K + k2);
        
        x1 = X_extended(:, k2);
        x2 = X_est(:, k1);
        x2 = align_to_reference_1D(x2, x1);
        
        plot(0:(L_optim-1), x1, 0:(L_optim-1), x2);
        
        if k1 == 1
            title(sprintf('True signal %d\n(blue)', k2));
        end
        if k2 == 1
            ylabel(sprintf('Estimated signal %d\n(orange)', k1));
        end
        
    end
end

savefig(gcf, sprintf('example_heterogeneous_%s_n_%d.fig', datestr(now(), 30), n));

%% Try to extract dominant subsignal of length L in each estimated signal

X_est_L = zeros(L, K);
for k = 1 : K
    for s = 0 : (L_optim - 1)
        x = circshift(X_est(:, k), s);
        x = x(1:L);
        if norm(x) > norm(X_est_L(:, k))
            X_est_L(:, k) = x;
        end
    end
end

figure(2);
clf;

for k1 = 1 : K
    for k2 = 1 : K
        subplot(K, K, (k1-1)*K + k2);
        
        x1 = X(:, k2);
        x2 = X_est_L(:, k1);
%         x2 = align_to_reference_1D(x2, x1);
        
        plot(0:(L-1), x1, 0:(L-1), x2);
        
        if k1 == 1
            title(sprintf('True signal %d\n(blue)', k2));
        end
        if k2 == 1
            ylabel(sprintf('Estimated signal %d\n(orange)', k1));
        end
        
    end
end

savefig(gcf, sprintf('example_heterogeneous_short_%s_n_%d.fig', datestr(now(), 30), n));

%% Reoptimize from the shortened estimator

gamma0 = gamma_est*(L/L_optim);

[X_est, gamma_est, problem2, stats2] = least_squares_1D_heterogeneous(moments, L, K, sigma_est, X_est_L, gamma0(:)); % check if the code fixes gamma to gamma0 or not


%% Redraw

fprintf('Estimated densities:\n');
disp(gamma_est');
fprintf('True densities:\n');
disp(m_actual*L/n);

figure(3);
clf;

for k1 = 1 : K
    for k2 = 1 : K
        subplot(K, K, (k1-1)*K + k2);
        
        x1 = X(:, k2);
        x2 = X_est(:, k1);
%         x2 = align_to_reference_1D(x2, x1);
        
        plot(0:(L-1), x1, 0:(L-1), x2);
        
        if k1 == 1
            title(sprintf('True signal %d\n(blue)', k2));
        end
        if k2 == 1
            ylabel(sprintf('Estimated signal %d\n(orange)', k1));
        end
        
    end
end

savefig(gcf, sprintf('example_heterogeneous_reoptimized_%s_n_%d.fig', datestr(now(), 30), n));


%%
permutations = perms(1:K);
for p = 1 : size(permutations, 1)
    P = permutations(p, :);
    fprintf('==\n');
    fprintf('Relative error subsignals of length L: %g\n', norm(X-X_est_L(:, P), 'fro') / norm(X, 'fro'));
    fprintf('Relative error after reoptimization:   %g\n', norm(X-X_est(:, P), 'fro') / norm(X, 'fro'));
end
fprintf('==\n');

save(sprintf('data_after_optim_example_heterogeneous_%s_n_%d.mat', datestr(now(), 30), n));
