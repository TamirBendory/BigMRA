clear; close all; clc;

%% Defining the problem

L = 5;
W = 2*L-1;
X = double(rgb2gray(imread('einstein_tongue_cropped.jpg')));
X = imresize(X, [L, L]);
xmin = min(X(:));
X = X - xmin;
xmax = max(X(:));
X = X / xmax;

sigma = 0.1;
m = 2000; % desired number of occurrences of the signal per micrograph
num_micrographs = 5; % number of micrographs
N = 2000; % size of single micrograph, in pixels on each side

%% generating data and computing moments
assert((W-1)/2 == round((W-1)/2), 'W assumed odd in this code.');
M1 = 0;
M2 = zeros(W);
m_eff = zeros(num_micrographs, 1);

for iter = 1:num_micrographs
    tic;
    [Y_clean, m_eff(iter)] = generate_clean_micrograph_2D(X, W, N, m);
    Y_obs = Y_clean + sigma*randn(N, N);
    fprintf('Micrograph %d generation: %.2g [s]\n', iter, toc());
    fprintf('   m_eff: %d\n', m_eff(iter));
    
    tic;
    M1_micrograph = sum(Y_obs(:));
    M2_micrograph = Second_moment_from_data(Y_obs, W, sigma);
    fprintf('Moment computation on micrograph %d: %.2g [s]\n', iter, toc());
    
    M1 = M1 + M1_micrograph;
    M2 = M2 + M2_micrograph;
    
end

% zero-padding the signal and compute moments for comparison
Xzp = zeros(W);
Xzp(1:L,1:L) = X;
PSX = abs(fft2(Xzp)).^2;
%MX1 = sum(X(:));

m_eff = sum(m_eff);
M1 = M1/m_eff;
M2 = M2/m_eff;
Px_est = real(fft2(M2));

PS_err = norm(Px_est(:) - PSX(:))/norm(PSX(:));
fprintf('PS estimation error = %.2g\n', PS_err);

%% estimating the signal using the RRR algorithm
tic
fprintf('RRR begins\n');
Xest = RRR(sqrt(Px_est),L);
Xest = Xest*sign(M1); % sign correction
fprintf('RRR %.2g [s]\n'),toc();

%%  plotting

figure;
imagesc([Xest, X]);

err = min(norm(Xest - X,'fro'), norm(rot90(Xest, 2) - X,'fro'))/norm(X(:));
fprintf('error = %.4g\n',err);
