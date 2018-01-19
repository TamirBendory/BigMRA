clear;
close all;
clc;

%% Defining the problem

% Load a grayscale image of size LxL and scale between 0 and 1.
L = 15;
W = 2*L-1;
X = double(rgb2gray(imread('einstein_tongue_cropped.jpg')));
X = imresize(X, [L, L]);
xmin = min(X(:));
X = X - xmin;
xmax = max(X(:));
X = X / xmax;

sigma = 0.0;
m = 2;
N = 100;

% if isempty(gcp('nocreate'))
%     parpool(2, 'IdleTimeout', 240);
% end

%% Generating data
tic;
[Y_clean, m_eff] = generate_clean_micrograph_2D(X, W, N, m);
Y_obs = Y_clean + sigma*randn(N, N);
fprintf('Gen data time: %.4g [s]\n', toc());
SNR = norm(Y_clean, 'fro')/norm(Y_obs-Y_clean, 'fro');
fprintf('SNR: %.4g\n', SNR);
fprintf('m_eff: %d\n', m_eff);

%% Pick which correlation coefficients de sample -- all for now
assert((W-1)/2 == round((W-1)/2), 'W assumed odd in this code.');

% Load lists of distinct moments of order 2 and 3 (if those were
% precomputed for that W; if not, call generate_list2_list3 with proper
% values of W (no need to recompute the ones that are already done: it
% takes a while.)
data = load(sprintf('lists_W_%d.mat', W));
assert(W == data.W);
list2 = data.list2;
list3 = data.list3;

% Can subsample if necessary
n3 = size(list3, 1);
nkeep = round(n3/W); % keep only nkeep elements
% % % nkeep = 0;
% nkeep = 10*W^2;
% nkeep = n3;
keep = randperm(n3, nkeep);
list3 = list3(keep, :);
n3 = size(list3, 1);


%%

tic;
[M1, M2, M3] = moments_from_data_no_debias_2D(Y_obs, list2, list3);
fprintf('Moment computation on micrograph: %.4g [s]\n', toc());

%%
X0 = [];
% X0 = X_zp; % cheat by giving true signal as initial guess
[X_est, problem, stats] = least_squares_2D(M1, M2, M3, W, sigma, N, L, m_eff, list2, list3, X0);

% Actually, BFGS seems to works nicely for this problem; keep it in mind.
% problem.linesearch = @(in1, in2) 2;
% X_est_2 = rlbfgs(problem);

%% Align X_est to X_zp (zero padded ground truth) for display
X_zp = [X zeros(L, size(X_est, 2)-L) ; zeros(size(X_est, 1)-L, size(X_est, 2))];
X_est = align_to_reference(X_est, X_zp);

imagesc([X_zp, X_est]); axis equal; axis tight; colormap gray;
savefig('latest.fig');