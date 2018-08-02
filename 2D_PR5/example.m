clear;
close all;
clc;

seed_rng = rng('shuffle');

if isempty(gcp('nocreate'))
    parpool(72,'IdleTimeout', 240);
end

%% Defining the problem

% Load a grayscale image of size LxL, removing the mean and scale between -1 and 1.
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));

L = 50;
W = 2*L-1;
X = imresize(X, [L, L]);
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
PSX = abs(fft2(X_zp)).^2;

% Pick a noise level
sigma = 3;
% Desired number of occurrences of the signal X in each micrograph
m_want = 1000;
%
micrograph_batch = 100; %512;
% Number of epocs
epocs = 1000;
% Number of micrographs to generate
n_micrographs = epocs*micrograph_batch;
% Each micrograph is square of size NxN
N = 4096;

err_PS = zeros(epocs,1);
err_rrr = zeros(epocs,1);

%% parameters for the RRR

th = 0.01;
%X_N = double(rgb2gray(imread('newton.jpg')));
%X_N = X_N/max(X_N(:));
%X_N = imresize(X_N, [L, L]);
%X_N = X_N - mean(X_N(:));
X_init = zeros(W);
%X_init(1:L,1:L) = X_N;
r = randn(L);
r = r - mean(r(:));
X_init(1:L,1:L) = r;

%% Generate the micrographs and collect their moments

%m_eff = zeros(n_micrographs, 1);
m = zeros(epocs,1);

for iter = 1:epocs
    
    % Generate micrographs
    Y_obs = zeros(N,N,micrograph_batch);
    tic;
    m_eff = zeros(micrograph_batch,1);
    for i = 1:micrograph_batch
        [Y_clean, m_eff(i)] = generate_clean_micrograph_2D(X, W, N, m_want);
        Y_obs(:,:,i) = Y_clean + sigma*randn(N);
    end
    fprintf('iter = %d,m  =  %d, time = %.4g [sec] \n', iter, sum(m_eff), toc);
    tic
    AC = ifft2(abs(fft2(Y_obs)).^2); %./m_eff;
    m(iter) = sum(m_eff);
    M2_micrograph = zeros(2*L-1,2*L-1,micrograph_batch);
    M2_micrograph(1:L,1:L,:) = AC(1:L,1:L,:);
    M2_micrograph(L+1:2*L-1,1:L,:) = AC(N-L+2:N,1:L,:);
    M2_micrograph(1:L,L+1:2*L-1,:) = AC(1:L,N-L+2:N,:);
    M2_micrograph(L+1:2*L-1,L+1:2*L-1,:) = AC(N-L+2:N,N-L+2:N,:);
    M2_micrograph(1,1,:) = M2_micrograph(1,1,:) - N^2*sigma^2;
    Normalize_factor = repmat(m_eff,[1,W,W]);
    Normalize_factor = permute(Normalize_factor,[2,3,1]);
    M2_micrograph = M2_micrograph./Normalize_factor;
    M2_micrograph = sum(M2_micrograph,3)./micrograph_batch;
    toc()
    if iter == 1
        M2 =  M2_micrograph;
    else
        M2 = M2 + M2_micrograph;
    end
    
    PS = fft2(M2/iter);
    err_PS(iter) = norm(PSX(:) - PS(:))/norm(PSX(:));
    fprintf('error PS = %.6g\n',err_PS(iter));
    max_iter = 200;
    [Xest_rrr, discrepancy_norm,err,err1, err2] = RRR_with_err(sqrt(PS),L,th,X_init,X,max_iter);
    [err_rrr(iter),min_iter] = min(err);
    fprintf('error RRR = %.4g\n',err_rrr(iter));
    save(strcat('Xest_',num2str(iter),'.mat'),'Xest_rrr');
    save('err_PS','err_PS');
    save('err_rrr','err_rrr'); 
end

   
