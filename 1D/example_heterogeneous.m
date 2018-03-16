clear;
close all;
clc;

%% Pick parameters and generate signals

% Pick K signals of length L and the size W of the separation window
K = 2;
L = 30;
W = 2*L-1;
X = randn(L, K);  % randn instead of rand here -- to be discussed

% Pick a noise level
sigma = .2;

% Desired number of occurrences of each signal X(:, k)
m_want = 10000*ones(K, 1);

% Length of micrograph
n = sum(m_want)*W*20;


%% Pick which correlation coefficients to sample

% Load lists of distinct moments of order 2 and 3
list2 = (0 : (L-1))';
% list2 = (1 : (L-1))';  % this one ignores biased terms

% Sampling strategy: see notes NB notebook 33
list3 = zeros(0, 2);
n3 = 0;
for k1 = 0 : (L-1)
    for k2 = 0 : -1 : -(L-1)
        if k1 - k2 <= L-1
            n3 = n3 + 1;
            list3(n3, :) = [k1, k2];
        end
    end
end

% Optionally, remove all moments that are affected by bias, so that it is
% no longer necessary to know sigma.
remove_biased_terms = true;
if remove_biased_terms
    list2(list2 == 0) = [];
    list3(list3(:, 1) == 0 | list3(:, 2) == 0, :) = [];
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

% Here we can choose the length of sought signals in optimization.
% Should be no less than L; typically set to W.
L_optim = 2*L-1;

X0 = randn(L_optim, K);
gamma0 = m_actual*L_optim/n; % give true gamma for now

[X_est, gamma_est, problem, stats] = least_squares_1D_heterogeneous(moments, L_optim, K, sigma, X0, gamma0(:)); % check if the code fixes gamma to gamma0 or not

% Now, should extract signals of length L out of the estimated signals of
% length L_optim and reoptimize (ideally).


%% Display

fprintf('Estimated densities:\n');
disp(gamma_est');
fprintf('True densities:\n');
disp(m_actual'*L_optim/n);

X_extended = [X ; zeros(L_optim-L, K)];

for k1 = 1 : K
    for k2 = 1 : K
        subplot(K, K, (k1-1)*K + k2);
        
        x1 = X_extended(:, k2);
        x2 = X_est(:, k1);
        x2 = align_to_reference_1D(x2, x1);
        
        plot(1:L_optim, x1, 1:L_optim, x2);
        
        if k1 == 1
            title(sprintf('True signal %d\n(blue)', k2));
        end
        if k2 == 1
            ylabel(sprintf('Estimated signal %d\n(orange)', k1));
        end
        
    end
end
