clear; close all; clc;
% script to generate progress figures for the paper

ind = [2,20,200];
%newton = rgb2gray(imread('newton.jpg'));
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
X = imresize(X, [50, 50]);

%figure(1); imagesc(newton); colormap gray; axis tight square off
%pdf_print_code(gcf, 'newton', 12)

for i = 1:length(ind)
str = strcat('images2\Xest_',num2str(ind(i)),'.mat');
load(str);
err1 = norm(Xest_rrr - X,'fro')/norm(X(:));
err2 = norm(rot90(Xest_rrr,2) - X,'fro')/norm(X(:));
if err2<err1
    Xest_rrr = rot90(Xest_rrr,2);
end
figure; imagesc(Xest_rrr); colormap gray; axis tight square off
str = strcat('reconstruction',num2str(ind(i)));
%pdf_print_code(gcf, str, 12)
end

%% progress

str = strcat('images2\err_rrr.mat');
load('err_rrr');
str = strcat('images2\err_PS.mat');
load('err_PS');
last_ind = max(find(err_rrr>0));
figure(11); 
subplot(121); loglog((1:last_ind)*500,err_PS(1:last_ind)); 
ylabel('PS estimation error');
xlabel('# micrographs')
xlim([500,500*last_ind])
axis square
grid on
subplot(122); 
loglog((1:last_ind)*500,err_rrr(1:last_ind)); 
ylabel('recovery error');
xlabel('# micrographs')
xlim([500,500*last_ind])
axis square
grid on

%pdf_print_code(gcf, 'Einstein_progress', 12)
