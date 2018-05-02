clear;
close all;
clc;

%% Pick parameters and generate signals

% Pick K signals of length L and the size W of the separation window
K = 3;
L = 21;
W = 2*L-1;
X = zeros(L, K);
X(:, 1) = [ones(ceil(L/2), 1) ; -ones(floor(L/2), 1)];
X(:, 2) = [linspace(0, 2, ceil(L/2))' ; linspace(2, 0, floor(L/2))'];
X(:, 3) = randn(L, 1);

% Pick a noise level
sigma = 3;

% Desired number of occurrences of each signal X(:, k)
m_want = [3e7 2e7 1e7];

% Length of micrograph
n = sum(m_want)*W*10;

fprintf('Micrograph length: %g\n\n\n', n);

%% Pick which correlation coefficients to sample

[list2, list3] = moment_selection(L, 'exclude biased');

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

L_optim = 2*L-1;
sigma_est = 0; % irrelevant if biased terms are excluded and if the weights internally do not depend on sigma

[X2, gamma2, X1, gamma1, X1_L, cost_X2] = heterogeneous_1D(moments, K, L, L_optim, sigma_est);

%%
permutations = perms(1:K);
for p = 1 : size(permutations, 1)
    P = permutations(p, :);
    fprintf('==\n');
    fprintf('Relative error subsignals of length L: %g\n', norm(X-X1_L(:, P), 'fro') / norm(X, 'fro'));
    fprintf('Relative error after reoptimization:   %g\n', norm(X-X2(:, P), 'fro') / norm(X, 'fro'));
end
fprintf('==\n');

% TODO: pick best P and apply it below

fprintf('Estimated densities:\n');
disp(gamma2');
fprintf('Estimated densities (before re optimization):\n');
disp(gamma1' * (L/L_optim));
fprintf('True densities:\n');
disp(m_actual*L/n);
