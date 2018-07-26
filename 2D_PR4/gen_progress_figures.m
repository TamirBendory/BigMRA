clear; close all; clc;
% script to generate progress figures for the paper

save_pdf = 0;
ind = [1,10,100, 480];
%newton = rgb2gray(imread('newton.jpg'));
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
X = imresize(X, [50, 50]);

%figure(1); imagesc(newton); colormap gray; axis tight square off
if save_pdf
pdf_print_code(gcf, 'newton', 12)
end

figure; 

for i = 1:length(ind)
str = strcat('images2\Xest_',num2str(ind(i)),'.mat');
load(str);
err1 = norm(Xest_rrr - X,'fro')/norm(X(:));
err2 = norm(rot90(Xest_rrr,2) - X,'fro')/norm(X(:));
if err2<err1
    Xest_rrr = rot90(Xest_rrr,2);
end
subplot(2,2,i); imagesc(Xest_rrr); colormap gray; axis tight square off
title(['# micrographs =', num2str(num2str(500*ind(i)))]);
if save_pdf
str = strcat('reconstruction',num2str(ind(i)));
pdf_print_code(gcf, str, 12)
end
end

%% progress

load('err_rrr');
load('err_PS');
last_ind = max(find(err_rrr>0));

figure(11); 
subplot(121); loglog((1:last_ind)*500,err_PS(1:last_ind)); 
ylabel('PS estimation error');
xlabel('# micrographs')
xlim([500,500*last_ind])
axis square
grid on

subplot(122); loglog((1:last_ind)*500,err_rrr(1:last_ind)); 
ylabel('recovery error');
xlabel('# micrographs')
xlim([500,500*last_ind])
axis square
grid on
if save_pdf
pdf_print_code(gcf, 'Einstein_progress', 12)
end