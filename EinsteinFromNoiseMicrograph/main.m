clear; clc; close all;

% Einstein's image
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
L = 50;
W = 2*L - 1;
X = imresize(X, [L, L]);
N = 4096;
m = 500;
num_micro = 1;
X_est = zeros(L,L,m);
sigma = 3;

for i = 1:num_micro
    [Y, placed] = generate_clean_micrograph_2D(X, W, N, m);
    Y = Y + sigma*randn(N);
    a = xcorr2(Y,X);
    [max_a, imax] = maxk((a(:)),placed);
    [xpeak, ypeak] = ind2sub(size(a),imax);
    for j = 1:placed
        X_est(:,:,j) = Y(xpeak(j)-L+1:xpeak(j),ypeak(j)-L+1:ypeak(j));
    end
end

Xest = mean(X_est,3);

figure; imagesc(Xest);
colormap gray; axis equal tight off square
%str = strcat('Einstein_from_noise_n',num2str(n));
%pdf_print_code(gcf, str , 12);

% figure;
% subplot(121); imagesc(X); colormap gray
% subplot(122); imagesc(Xest); colormap gray
