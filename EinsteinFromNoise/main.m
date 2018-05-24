clear; clc; close all; 

% Einstein's image
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
L = 50; 
X = imresize(X, [L, L]);

% number of images
n = 1000000; 
%N = circshift(X,[10,20]);
Xest = zeros(L,L,10);
for j = 1:10
 j
N_aligned = zeros(L,L,n/10);
N = 0.4578*randn(L,L,n/10);
for i = 1:n/10
    a = xcorr2(X,N(:,:,i));
    [max_a, imax] = max((a(:)));
    [xpeak, ypeak] = ind2sub(size(a),imax(1));
   % N_aligned(:,:,i) = 
    N_aligned(:,:,i) = circshift(N(:,:,i),[xpeak-L,ypeak-L]);
end
Xest(:,:,j) = mean(N_aligned,3);
end

Xest1 = mean(Xest,3);

figure; imagesc(Xest1);
colormap gray; axis equal tight off square
str = strcat('Einstein_from_noise_n',num2str(n));
pdf_print_code(gcf, str , 12);

% figure; 
% subplot(121); imagesc(X); colormap gray
% subplot(122); imagesc(Xest); colormap gray
