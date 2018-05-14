clear;
close all;
clc;

seed_rng = rng('shuffle');

%if isempty(gcp('nocreate'))
%    parpool(20,'IdleTimeout', 240);
%end

%% Defining the problem

% Load a grayscale image of size LxL, removing the mean and scale between -1 and 1.
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
L = 50;
W = 2*L-1;
X = imresize(X, [L, L]);
% Pick a noise level
sigma = 3;
% Desired number of occurrences of the signal X in each micrograph
m_want = 1000; 
% Number of micrographs to generate
n_micrographs = 1; %1e4;
% Each micrograph is square of size NxN
N = 4000; %1000;

%% Pick which correlation coefficients to sample

top = (W-1)/2;
list2 = zeros((W^2+1)/2, 2);
k = 0;
for k1 = 0 : top
    k = k + 1;
    list2(k, :) = [k1, 0];
end
for k1 = (-top) : top
    for k2 = 1 : top
        k = k + 1;
        list2(k, :) = [k1, k2];
    end
end

list3 = [];
n3 = size(list3, 1);
%% Generate the micrographs and collect their moments

m_eff = zeros(n_micrographs, 1);

for iter = 1:n_micrographs

% Generate micrograph
    tic;
    [Y_clean, m_eff(iter)] = generate_clean_micrograph_2D(X, W, N, m_want);
    Y_obs = Y_clean + sigma*randn(N);    
    %Y_obs = imgaussfilt(Y_obs,4);
    fprintf('iter = %d, m_eff =  %d \n', iter, m_eff(iter));
    [M1_micrograph, M2_micrograph, M3_micrograph] = moments_from_data_no_debias_2D_v5(Y_obs, list2, list3);    
    fprintf('Iteration timing: %.2g [s]\n', toc());
    
    if iter == 1
        M1 =  M1_micrograph;
        M2 =  M2_micrograph;
        M3 =  M3_micrograph;
    else
        M1 = M1 + M1_micrograph;
        M2 = M2 + M2_micrograph;
        M3 = M3 + M3_micrograph;
    end
    
    if mod(iter,100) == 0
        % removing big and unnecessary variables before saving  
        clear 'Y_clean' 'Y_obs' 'list2' 'M2_micrograph'
        save(strcat('data_exp_',num2str(iter),'_latte.mat'));
    end
    
end
