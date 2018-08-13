clear; close all; clc;
% script to generate progress figures for the paper

save_pdf = 0;
ind = [1,10,100,600];
%newton = rgb2gray(imread('newton.jpg'));
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
X = imresize(X, [50, 50]);

figure; 

for i = 1:length(ind)
str = strcat('images\Xest_',num2str(ind(i)),'.mat');
load(str);
Xest_rrr = Xest_rrr(1:size(X,1),1:size(X,2));
err1 = norm(Xest_rrr - X,'fro')/norm(X(:));
err2 = norm(rot90(Xest_rrr,2) - X,'fro')/norm(X(:));
if err2<err1
    Xest_rrr = rot90(Xest_rrr,2);
end
subplot(2,2,i); imagesc(Xest_rrr); colormap gray; axis tight square off
title(['# micrographs = 10^', num2str(log10(100*ind(i)))]);
% if save_pdf
% str = strcat('reconstruction',num2str(ind(i)));
% pdf_print_code(gcf, str, 12)
% end
end

if save_pdf
str = strcat('Einstien_progress_examples');
pdf_print_code(gcf, str, 12)
end

%% progress

% load('err_rrr');
% load('err_PS');
load('images\err_rrr');
load('images\err_PS');

last_ind = max(find(err_rrr(:,end)>0));

figure(11); 
subplot(121); loglog((1:last_ind)*100,err_PS(1:last_ind),'.b'); 
ylabel('PS estimation error');
xlabel('# micrographs')
xlim([100,100*last_ind])
axis square
grid on

subplot(122); loglog((1:last_ind)*100,err_rrr(1:last_ind),'.b'); 
ylabel('recovery error');
xlabel('# micrographs')
xlim([100,100*last_ind])
%xlim([10^2,10^5])
ylim([.2,.7])

axis square
grid on
if save_pdf
pdf_print_code(gcf, 'Einstein_progress', 12)
end