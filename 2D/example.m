clear;
close all;
clc;

%% Defining the problem

% Load a grayscale image of size LxL and scale between 0 and 1.
L = 3;
W = 2*L-1;
X = double(rgb2gray(imread('einstein_tongue_cropped.jpg')));
X = imresize(X, [L, L]);
xmin = min(X(:));
X = X - xmin;
xmax = max(X(:));
X = X / xmax;

sigma = 0.1;
m = 100;
N = 800;

if isempty(gcp('nocreate'))
    parpool(2, 'IdleTimeout', 240);
end

%% Generating data
tic;
[Y_obs, Y_clean, ind, class] = gen_data2D(X, N, m, sigma, W);
fprintf('Gen data time: %.2g [s]\n', toc());
snr = norm(Y_clean(:))/norm(Y_obs(:)-Y_clean(:));
m_eff = size(ind, 1);
fprintf('SNR: %.2g\n', snr);
fprintf('m_eff: %d\n', m_eff);

%% Pick which correlation coefficients de sample -- all for now
list2 = zeros(W^2, 2);
k = 0;
range = (-(W-1)/2) : ((W-1)/2);
for k1 = range
    for k2 = range
        k = k + 1;
        list2(k, :) = [k1, k2];
    end
end

list3 = zeros(W^4, 4);
k = 0;
for k1 = range
    for k2 = range
        for l1 = range
            for l2 = range
                k = k + 1;
                list3(k, :) = [k1, k2, l1, l2];
            end
        end
    end
end


%%
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];

[M1, M2, M3] = moments_from_data_no_debias_2D(Y_obs, list2, list3);

X0 = [];
% X0 = X_zp; % cheat by giving true signal as initial guess
[X_est, problem] = least_squares_2D(M1, M2, M3, W, sigma, N, L, m_eff, list2, list3, []);


% Align X_est to X_zp for display
X_est = align_to_reference(X_est, X_zp);

imagesc([X_zp, X_est]); axis equal;