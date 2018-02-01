clear;
close all;
clc;

%% Defining the problem

% Pick a signal of length L and the size W of the separation window
L = 6;
W = 2*L-1;
x = rand(L, 1); % rand, not randn (to be discussed)

% Pick a noise level
sigma = .1;

% Desired number of occurrences of the signal x in each micrograph
m_want = 100;
% Number of micrographs to generate
n_micrographs = 1;
% Each micrograph has length N
N = m_want*W*10;

%% Pick which correlation coefficients to sample

assert((W-1)/2 == round((W-1)/2), 'W assumed odd in this code.');

% Load lists of distinct moments of order 2 and 3
list2 = (0 : (W-1)/2)';
range = (-(W-1)/2 : (W-1)/2)';
[p, q] = meshgrid(range);
list3 = [p(:), q(:)]; % there is redudance in here.

% Can subsample if necessary:  keep only nkeep triple correlation coeffs.
% n3 = size(list3, 1);
% nkeep = round(n3/W);
% keep = randperm(n3, nkeep);
% list3 = list3(keep, :);

n3 = size(list3, 1);

%% Generate the micrographs and collect their moments

m_eff = zeros(n_micrographs, 1);

for iter = 1:n_micrographs
    
    % Generate micrograph
    tic;
    [y_clean, m_eff(iter)] = generate_clean_micrograph_1D(x, W, N, m_want);
    y_obs = y_clean + sigma*randn(N, 1);
    fprintf('Time to generate micrograph %3d: %.2g [s]\n', iter, toc());
    SNR = norm(y_clean, 'fro')/norm(y_obs-y_clean, 'fro');
    fprintf('   SNR: %.2g\n', SNR);
    fprintf('   m_eff: %d\n', m_eff(iter));
    
    % Compute the moments and aggregate
    tic;
    [M1_micrograph, M2_micrograph, M3_micrograph] = moments_from_data_no_debias_1D(y_obs, list2, list3);
    fprintf('   Moment computation: %.2g [s]\n', toc());
    
    if iter == 1
        M1 =  M1_micrograph;
        M2 =  M2_micrograph;
        M3 =  M3_micrograph;
    else
        M1 = M1 + M1_micrograph;
        M2 = M2 + M2_micrograph;
        M3 = M3 + M3_micrograph;
    end
    
end

% Total number of occurrences of X in the micrographs
m_total = sum(m_eff);

%% Optimization

% TODO: think about the equivalent micrograph size when we have multiple
% micrographs. This doesn't sound right to me. For sigma = 0, it's
% irrelevant because N only intervenes in the bias computations.
N_eff = N*n_micrographs;

x0 = [];
[x_est, problem] = least_squares_1D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, x0);

%% Display

% Align x_est to x_zp (zero padded ground truth) for display
x_zp = [x ; zeros(W-L, 1)];
x_est = align_to_reference_1D(x_est, x_zp);
% X_est_aligned = align_by_energy(X_est,L);

% err = norm(X(:) -  X_est_aligned(:))/norm(X(:));
err = norm(x_zp(:) -  x_est(:))/norm(x(:));
fprintf('error = %.4g\n',err);
% figure(1); imagesc([X, X_est_aligned]); axis equal;
figure(1);
plot(1:W, x_zp, '.-', 1:W, x_est, 'o-');
legend('Truth', 'Estimated');
