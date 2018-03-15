clear;
close all;
clc;

%% Pick parameters and generate signals

% Pick K signals of length L and the size W of the separation window
K = 1;
L = 6;
W = 2*L-1;
X = rand(L, K);  % randn instead of rand here -- to be discussed

% Pick a noise level
sigma = 0;

% Desired number of occurrences of each signal X(:, k)
m_want = 100*ones(K, 1);

% Length of micrograph
n = sum(m_want)*W*10;


%% Pick which correlation coefficients to sample

% Load lists of distinct moments of order 2 and 3
list2 = (0 : (L-1))';
% list2 = (1 : (L-1))';  % this one ignores biased terms

list3 = zeros(0, 2);
n3 = 0;
for k1 = 0 : (L-1)
    for k2 = 0 : (L-1)
        if k1 + k2 >= L
            continue;
        end
        n3 = n3 + 1;
        list3(n3, :) = [k1, k2];
    end
end

%% Generate the micrograph

tic;
[y_clean, m_actual] = generate_clean_micrograph_1D_heterogeneous(X, W, n, m_want);
y_obs = y_clean + sigma*randn(n, 1);
fprintf('Time to generate micrograph: %.2g [s]\n', toc());
SNR = norm(y_clean, 'fro')/norm(y_obs-y_clean, 'fro');
fprintf('   SNR: %.2g\n', SNR);
fprintf('   m_actual/m_want: ');
fprintf(' %.2g', m_actual./m_want);
fprintf('\n');

% Total number of occurrences of signals in the micrograph
m_total = sum(m_actual);

%% Collect the moments
tic;
[M1, M2, M3] = moments_from_data_no_debias_1D(y_obs, list2, list3);
fprintf('   Moment computation: %.2g [s]\n', toc());

moments.M1 = M1 / n;  %%%% !!!!! We normalize by n here; ideally, we should normalize by n in moments_from_data_no_debias_1D, but I'll only do the change when everything is under control
moments.M2 = M2 / n;
moments.M3 = M3 / n;
moments.list2 = list2;
moments.list3 = list3;

%% Optimization

[X_est_W, gamma_est, problem, stats] = least_squares_1D_heterogeneous(moments, W, K, sigma);

% Now, extract signals of length L out of the estimated signals of length W


%% Display
plot(X_est_W);
hold all;
plot(X);

gamma_est
m_actual*W/n
