clear;
close all;
clc;

%% Pick parameters and generate signal

% Pick one signal of length L and the size W of the separation window.
L = 21;
W = 2*L-1;
X = randn(L, 1);

% Pick a noise level.
sigma = 2;

% Desired length of micrograph.
n = 1e9;

% Desired number of occurrences of the signal.
m_want = floor(n/(5*W));

fprintf('Micrograph length: %g\n\n\n', n);

%% Pick which correlation coefficients to sample

[list2, list3] = moment_selection(L, 'exclude biased');

%% Generate the micrograph

[y_clean, m_actual] = generate_clean_micrograph_1D_heterogeneous(X, W, n, m_want);
y_obs = y_clean + sigma*randn(n, 1);
gamma = m_actual*L/n;
SNR = norm(y_clean, 'fro')/norm(y_obs-y_clean, 'fro');

%% The grand experiment starts here

% Select sizes of sub-micrographs to consider.
ns = unique(round(logspace(4, log10(n), 9)));

% How many times do we optimize from a different random initial guess?
n_init_optim = 3;

for iter = 1 : length(ns)

    % Collect the moments for the first bit of the micrograph.
    nn = ns(iter);
    [M1, M2, M3] = moments_from_data_no_debias_1D_batch( ...
                                           y_obs(1:nn), list2, list3, 1e8);
	%! We normalize by nn here.
    moments.M1 = M1 / nn;
    moments.M2 = M2 / nn;
    moments.M3 = M3 / nn;
    moments.list2 = list2;
    moments.list3 = list3;

    % Parameters and initializations for optimization:
    % empty inputs mean default values are picked.
    L_optim = 2*L-1;
    sigma_est = [];
    X0 = [];
    gamma0 = []; % True value is: m_actual*L_optim/n
    
    % Run the optimization from different random initializations n_repeat
    % times, and keep the best result according to the cost value of X2.
    result = struct();
    result.cost_X2 = inf;
    for repeat = 1 : n_init_optim
        [X2, gamma2, X1, gamma1, X1_L, cost_X2] = heterogeneous_1D( ...
                            moments, 1, L, L_optim, sigma_est, X0, gamma0);
        if cost_X2 < result.cost_X2
            result.X1 = X1;
            result.X1_L = X1_L;
            result.X2 = X2;
            result.gamma1 = gamma1;
            result.gamma2 = gamma2;
            result.RMSE1 = norm(X1_L-X)/norm(X);
            result.RMSE2 = norm(X2-X)/norm(X);
            result.cost_X2 = cost_X2;
        end
    end

    % Save in a structure array.
    results(iter) = result; %#ok<SAGROW>
    
end

clear y_clean y_obs;

ID = randi(1000000);

save(sprintf('progressive_n%d_%d.mat', n, ID));

%%
figure(1);

for iter = 1 : length(ns)
    result = results(iter);
    subplot(3, 3, iter);
    title(sprintf('n = %d', nn));
    T = 0:(L-1);
    plot(T, X, T, result.X1_L, T, result.X2);
%     hleg = legend(sprintf('Ground truth (%.2g)', result.gamma), ...
%                   sprintf('First estimate (%.2g)', result.gamma1*L/L_optim), ...
%                   sprintf('Final estimate (%.2g)', result.gamma2));
%     set(hleg, 'Location', 'northoutside');
%     set(hleg, 'Orientation', 'horizontal');
end

savefig(gcf, sprintf('progressive_n%d_%d.fig', n, ID));

%%
figure(2);

loglog(ns, [results.RMSE1], '.-', ns, [results.RMSE2], '.-');
hleg = legend('First estimate', 'Final estimate');
set(hleg, 'Location', 'northoutside');
set(hleg, 'Orientation', 'horizontal');
title('Root mean squared error');

savefig(gcf, sprintf('progressive_RMSE_n%d_%d.fig', n, ID));
