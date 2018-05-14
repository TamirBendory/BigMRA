clear;
%close all;
clc;
rng('shuffle')

%% Defining the problem

L = 52;

% Load a grayscale image of size LxL and scale between 0 and 1.
X_origin = rgb2gray(imread('Einstein5_small.jpg'));
X_origin = double(imresize(X_origin,[L,L]));
X = X_origin;
%X = X/max(X(:));
%X = X - imgaussfilt(X,5);
X = X/max(X(:));
X = X - mean(X(:));
if 1
figure; 
subplot(221); imagesc(X_origin); colormap gray
subplot(222); imagesc(X); colormap gray
subplot(223); imagesc(fftshift(abs(fft2(X_origin)))); colormap gray
subplot(224); imagesc(fftshift(abs(fft2(X)))); colormap gray
end
W = 2*L-1;
% Pick a noise level
sigma = 2;
% Desired number of occurrences of the signal X in each micrograph
m_want = 2; 
% Number of micrographs to generate
n_micrographs = 1;
% Each micrograph is square of size NxN
N = 400; %1000;
% LP filtering flag
LP_flag = 1;
% LP filtering std
sigma_gauss = 4;
% Over-crossing factor for noisy micrgraph
Over_cross_factor = 2;
%% Generate the micrographs and collect their moments

[Y_clean, m_eff] = generate_clean_micrograph_2D(X, W, N, m_want);
noise = sigma*randn(N);
Y_obs = Y_clean + noise;
fprintf('SNR = %.4g\n',norm(Y_clean(:))^2/norm(noise(:))^2);
if LP_flag
Y_obs = imgaussfilt(Y_obs,sigma_gauss);
end
%figure(1); imagesc(X); colormap gray; title('target image');
%figure(2); imagesc(Y_clean); colormap gray; title('microgrpah');

% noisy XC
X_xc = X; %+sigma*randn(L);
%if LP_flag
%X_xc = imgaussfilt(X_xc,sigma_gauss);
%end
Xc_noisy = xcorr2(X_xc,Y_obs);
Xc_noisy = rot90(Xc_noisy(L:end,L:end)/max(Xc_noisy(:)),2);
[max_valn max_indn]= maxk(Xc_noisy(:),Over_cross_factor*m_eff);
[ind2n, ind1n] = ind2sub(N,max_indn);
%Max_map_noisy = insertMarker(Y_obs,[ind1n ind2n],'x','color','green','size',14);
Max_map_noisy = Y_obs;
% clean XC
Xc_clean = xcorr2(X,Y_clean);
Xc_clean = rot90(Xc_clean(L:end,L:end)/max(Xc_clean(:)),2);
[max_val max_ind]= maxk(Xc_clean(:),m_eff);
[ind2, ind1] = ind2sub(N,max_ind);
color = '{';
for i = 1:m_eff
    color = strcat(color,'''red'',');
end
for i = 1:Over_cross_factor*m_eff
    if i == Over_cross_factor*m_eff
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
