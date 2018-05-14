clear;
close all;
clc;

rng('shuffle')

if isempty(gcp('nocreate'))
    parpool(32,'IdleTimeout', 240);
end

%% Defining the problem

% Load a grayscale image of size LxL and scale between 0 and 1.
X = double(rgb2gray(imread('einstein_tongue_free_backround.jpg')));
X = X/max(X(:));
L = 100;
W = 2*L-1;
X = imresize(X, [L, L]);
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
% Pick a noise level
sigma = 3;

% Desired number of occurrences of the signal X in each micrograph
m_want = 1000; %3000;
% Number of micrographs to generate
n_micrographs = 1e6;
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
    %    Y_obs = Y_clean + sigma*randn(N);
    Y_obs = Y_clean + sigma*randn(N);
    
    fprintf('iter = %d\n', iter);
    fprintf('m_eff: %d\n', m_eff(iter));
    
    % Compute the moments and aggregate
 %   tic;
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
        save('data_exp_chai.mat');
    end
    
end
