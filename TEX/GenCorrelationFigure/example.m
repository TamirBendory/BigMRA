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
m_want = 4; %1000;
% Number of micrographs to generate
%n_micrographs = 1; %1e4;
% Each micrograph is square of size NxN
N = 250; %4000; %1000;

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
% Generate micrograph
%[Y_clean, m_eff] = generate_clean_micrograph_2D(X, W, N, m_want);
Y_clean = zeros(N); 
ind1 = [20,25];
ind2 = [180, 35];
ind3 = [35, 145];
ind4 = [150, 160];

Y_clean (ind1(1):ind1(1)+L-1,ind1(2):ind1(2)+L-1) = X;
Y_clean (ind2(1):ind2(1)+L-1,ind2(2):ind2(2)+L-1) = X;
Y_clean (ind3(1):ind3(1)+L-1,ind3(2):ind3(2)+L-1) = X;
Y_clean (ind4(1):ind4(1)+L-1,ind4(2):ind4(2)+L-1) = X;

Y_obs = Y_clean + sigma*randn(N);

figure; imagesc(Y_clean); colormap gray; axis equal tight off
filename = 'micrograph_Einstein_example_clean.pdf';
pdf_print_code(gcf, filename, 11)
figure; imagesc(Y_clean + 0.5*randn(N)); colormap gray; axis equal tight off
filename = 'micrograph_Einstein_example_s05.pdf';
pdf_print_code(gcf, filename, 11)
figure; imagesc(Y_clean + 3*randn(N)); colormap gray; axis equal tight off
filename = 'micrograph_Einstein_example_s3.pdf';
pdf_print_code(gcf, filename, 11)
