clear; clc; close all; 

% Einstein's image
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
L = 50; 
X = imresize(X, [L, L]);

% number of images
n = 100000; 
N = 0.4578*randn(L,L,n);
%N = circshift(X,[10,20]);
N_aligned = zeros(L,L,n);

for i = 1:n
    a = xcorr2(X,N(:,:,i));
    [max_a, imax] = max((a(:)));
    [xpeak, ypeak] = ind2sub(size(a),imax(1));
   % N_aligned(:,:,i) = 
    N_aligned(:,:,i) = circshift(N(:,:,i),[xpeak-L,ypeak-L]);
end

Xest = mean(N_aligned,3);

figure; imagesc(Xest);
colormap gray; axis equal tight off square
str = strcat('Einstein_from_noise_n',num2str(n));
pdf_print_code(gcf, str , 12);

% figure; 
% subplot(121); imagesc(X); colormap gray
% subplot(122); imagesc(Xest); colormap gray
