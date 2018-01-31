clear;
close all;
clc;

%% Defining the problem

% Load a grayscale image of size LxL and scale between 0 and 1.
L = 4;
W = 2*L-1;
X = double(rgb2gray(imread('einstein_tongue_cropped.jpg')));
X = imresize(X, [L, L]);
xmin = min(X(:));
X = X - xmin;
xmax = max(X(:));
X = X / xmax;

% Pick a noise level
sigma = .01;

% Desired number of occurrences of the signal X in each micrograph
m_want = 3000;
% Number of micrographs to generate
n_micrographs = 1;
% Each micrograph is square of size NxN
N = 2000;

%% Pick which correlation coefficients to sample

assert((W-1)/2 == round((W-1)/2), 'W assumed odd in this code.');

% Load lists of distinct moments of order 2 and 3 (if those were
% precomputed for that W; if not, call generate_list2_list3 with proper
% values of W (no need to recompute the ones that are already done: it
% takes a while.)
data = load(sprintf('lists_W_%d.mat', W));
assert(W == data.W);
list2 = data.list2;
list3 = data.list3;

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
    [Y_clean, m_eff(iter)] = generate_clean_micrograph_2D(X, W, N, m_want);
    Y_obs = Y_clean + sigma*randn(N);
    fprintf('Time to generate micrograph %3d: %.2g [s]\n', iter, toc());
    SNR = norm(Y_clean, 'fro')/norm(Y_obs-Y_clean, 'fro');
    fprintf('   SNR: %.2g\n', SNR);
    fprintf('   m_eff: %d\n', m_eff(iter));
    
    % Compute the moments and aggregate
    tic;
    [M1_micrograph, M2_micrograph, M3_micrograph] = moments_from_data_no_debias_2D(Y_obs, list2, list3);
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
N_eff = N*sqrt(n_micrographs);

X0 = [];
[X_est, problem] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0);

%% Display

% Align X_est to X_zp (zero padded ground truth) for display
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
X_est = align_to_reference(X_est, X_zp);
% X_est_aligned = align_by_energy(X_est,L);

% err = norm(X(:) -  X_est_aligned(:))/norm(X(:));
err = norm(X_zp(:) -  X_est(:))/norm(X(:));
fprintf('error = %.4g\n',err);
% figure(1); imagesc([X, X_est_aligned]); axis equal;
figure(1);
imagesc([X_zp, X_est]);
axis equal;
