clear;
close all;
clc;

rng('shuffle')

%% Defining the problem

L = 52;

% Load a grayscale image of size LxL and scale between 0 and 1.
X = rgb2gray(imread('Einstein5_small.jpg'));
X = double(imresize(X,[L,L]));
X = X/max(X(:));
Xf = fft2(X);
W = 2*L-1;
% Pick a noise level
sigma = 0.005; %9;
% Desired number of occurrences of the signal X in each micrograph
m_want = 3; %3000;
% Number of micrographs to generate
n_micrographs = 1;
% Each micrograph is square of size NxN
N = 500; %1000;

%% Generate the micrographs and collect their moments

[Y_clean, m_eff] = generate_clean_micrograph_2D(X, W, N, m_want);
Yf = fft2(Y_clean);
n = randn(N);
noise = abs(Yf).*sigma.*fft2(n);
Y_obs_f = Yf+noise;
Y_obs = real(ifft2(Y_obs_f));
fprintf('SNR = %.4g\n',norm(Yf(:))^2/norm(noise(:))^2);
%Y_obs = imgaussfilt(Y_obs,2);

figure(1); imagesc(X); colormap gray; title('target image');
figure(2); imagesc(Y_clean); colormap gray; title('microgrpah');

% noisy XC
Xc_noisy = xcorr2(X+sigma*randn(L),Y_obs);
Xc_noisy = rot90(Xc_noisy(L:end,L:end)/max(Xc_noisy(:)),2);
[max_valn max_indn]= maxk(Xc_noisy(:),m_eff);
[ind2n, ind1n] = ind2sub(N,max_indn);
Max_map_noisy = insertMarker(Y_obs,[ind1n ind2n],'x','color','green','size',14);

% clean XC
Xc_clean = xcorr2(X,Y_clean);
Xc_clean = rot90(Xc_clean(L:end,L:end)/max(Xc_clean(:)),2);
[max_val max_ind]= maxk(Xc_clean(:),m_eff);
[ind2, ind1] = ind2sub(N,max_ind);
color = '{';
for i = 1:m_eff
    color = strcat(color,'''red'',');
end
for i = 1:m_eff
    if i == m_eff
    color = strcat(color,'''green''');
    else
        color = strcat(color,'''green'',');    
    end
end
color = strcat(color,'}');

Max_map_clean = insertMarker(Y_clean,[[ind1; ind1n] [ind2; ind2n]],'color',eval(color),'size',14);

figure(3);  
subplot(131); imshow(Max_map_clean); title('Clean')
subplot(132); imshow(Max_map_noisy); title('Noisy')
subplot(133); imshow(Xc_noisy); title('XC Noisy')
