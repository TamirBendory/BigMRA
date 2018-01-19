clear;
close all;
clc;

%% Defining the problem

% Load a grayscale image of size LxL and scale between 0 and 1.
% <<<<<<< HEAD
L = 21;
% =======
% L = 15;
% >>>>>>> 465debbc500803ffe1d3ce4e9add830d9f1f8595
W = 2*L-1;
X = double(rgb2gray(imread('einstein_tongue_cropped.jpg')));
X = imresize(X, [L, L]);
xmin = min(X(:));
X = X - xmin;
xmax = max(X(:));
X = X / xmax;
%X = X - mean(X(:));

%<<<<<<< HEAD
sigma = .5;
m = 200; % per micrograph
Num_micrographs = 100;
N = 1200; % size of single micrograph

% if isempty(gcp('nocreate'))
%     parpool(2, 'IdleTimeout', 480);
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
[list2, list3] = list_distinct_moments_2D(W,0);
%
% data = load(sprintf('lists_W_%d.mat', W));
% assert(W == data.W);
% list2 = data.list2;
% list3 = data.list3;

% Can subsample if necessary
n3 = size(list3, 1);
nkeep = round(n3/(W^2)); % keep only nkeep elements
keep = randperm(n3, nkeep);
list3 = list3(keep, :);
n3 = size(list3, 1);


m_eff = zeros(Num_micrographs,1);

M1 = 0;
M2 = zeros(size(list2,1),1);
M3 = zeros(size(list3,1),1);

for i = 1:Num_micrographs
    tic;
    [Y_clean, m_eff(i)] = generate_clean_micrograph_2D(X, W, N, m);
    Y_obs = Y_clean + sigma*randn(N, N);
    fprintf('Gen  micrograph %d time: %.2g [s]\n', i,toc());
    SNR = norm(Y_clean, 'fro')/norm(Y_obs-Y_clean, 'fro');
    fprintf('SNR: %.2g\n', SNR);
    fprintf('m_eff: %d\n', m_eff(i));
    
    %%
    
    tic;
    [M1_micrograph, M2_micrograph, M3_micrograph] = moments_from_data_no_debias_2D(Y_obs, list2, list3);
    fprintf('Moment computation on micrograph: %.2g [s]\n', toc());
    
    M1 = M1 + M1_micrograph;
    M2 = M2 + M2_micrograph;
    M3 = M3 + M3_micrograph;
    
end

%  M1 = M1/Num_micrograph;
%  M2 = M2/Num_micrograph;
%  M3 = M3/Num_micrograph;

m_eff = sum(m_eff);
%=======
% tic;
% [M1, M2, M3] = moments_from_data_no_debias_2D(Y_obs, list2, list3);
% fprintf('Moment computation on micrograph: %.4g [s]\n', toc());

%X0 = [];
% X0 = X_zp; % cheat by giving true signal as initial guess
[X_est, problem] = least_squares_2D(M1, M2, M3, W, sigma, round(N*sqrt(Num_micrographs)), L, m_eff, list2, list3, []);

% Actually, BFGS seems to works nicely for this problem; keep it in mind.
% problem.linesearch = @(in1, in2) 2;
% X_est_2 = rlbfgs(problem);

% Align X_est to X_zp (zero padded ground truth) for display
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
%X_est = align_to_reference(X_est, X_zp);
X_est_aligned = align_by_energy_xcorr(X_est,L);

%<<<<<<< HEAD
err   = norm(X(:) -  X_est_aligned(:))/norm(X(:));
%err   = norm(X_zp(:) -  X_est(:))/norm(X(:));
fprintf('error = %.4g\n',err);
figure(1); imagesc([X, X_est_aligned]); axis equal;
%=======
%imagesc([X_zp, X_est]); axis equal;
%savefig('latest.fig');
%>>>>>>> 465debbc500803ffe1d3ce4e9add830d9f1f8595
