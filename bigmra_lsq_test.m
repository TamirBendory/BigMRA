clear;
close all;
clc;

%% Problem parameters
L = 11;             % length of signal
sigma = 1;          % noise level
W = 2*L;            % window length
Nfactor = 4;        % Sparsity factor
K = 1;              % Number of distinct signals (heterogeneity)
m = 1e3*ones(K, 1); % Number of repetitions of each signal

%% Generate signals
x_true = randn(L, K);

%% Generate data
N = W*Nfactor*sum(m);
[y, yc, ind, class] = gen_data(x_true, N, m, sigma, W);
m_eff = zeros(K,1);
for k = 1:K
    m_eff(k) = sum(class == k);
end

%%  Process the data: compute moments

window_shift = 1;   % By how much does the window slide at each step?
y_mat = gen_data_mtx(y, W, window_shift);

if isempty(gcp('nocreate'))
    parpool(2, 'IdleTimeout', 240);
end

% computing the empirical invariants of the data
[M1, M2, M3] = invariants_from_data_no_debias(y_mat);

%% Run the least squares

x_est = bigmra_lsq(M1, M2, M3, W, sigma, N, L, m_eff);


%% Visualization
if K == 1
    plot(1:L, x_true, 'o-', 1:L, x_est, '.-');
    title(sprintf('Relative error: %g', norm(x_true - x_est, 'fro')/norm(x_true, 'fro')));
end