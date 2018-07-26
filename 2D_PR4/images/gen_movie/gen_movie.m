%%% generating a movie
clear;
close all;
clc;

% saving to video
v = VideoWriter('Einstein.avi');
v.FrameRate = 10;
open(v);
fig = figure('position',[100 100 500 500]);

% data parameters
max_ind = 842;
axisX = 512*(1:max_ind);
L = 50;
W = 2*L-1;
load('err_PS');
load('err_rrr');
err_PS = err_PS(1:max_ind);
err_rrr = err_rrr(1:max_ind);
err_rrr (140) = err_rrr (139);

% Load a grayscale image of size LxL, removing the mean and scale between -1 and 1.
X = double(rgb2gray(imread('Einstein5_small.jpg')));
X = X - mean(X(:));
X = X/max(abs(X(:)));
L = 50; 
X = imresize(X, [L, L]);

%% First frame

X_N = double(rgb2gray(imread('newton.jpg')));
X_N = X_N/max(X_N(:));
X_N = imresize(X_N, [L, L]);
X_N = X_N - mean(X_N(:));
X_init = zeros(W);
X_init(1:L,1:L) = X_N;
subplot(2,2,[1,2]); imagesc(X_N);
colormap gray; axis equal tight off square
shg;
Frame(1) = getframe(fig);
writeVideo(v,Frame(1));


%% movie
for i = 1:max_ind
    
    % loading data
    str = strcat('../Xest_',num2str(i));
    load(str);
    % would be removed in the future
    err1 = norm(Xest_rrr - X,'fro');
    err2 = norm(rot90(Xest_rrr,2) - X,'fro');
    if err2<err1
        Xest_rrr = rot90(Xest_rrr,2);
    end
    subplot(2,2,[1,2]); imagesc(Xest_rrr);
    colormap gray; axis equal tight off square
    
    
    subplot(2,2,3); loglog(axisX,[err_PS(1:i);zeros(max_ind-i,1)],'b');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('PS estimation error')
    xlim([512, max_ind*512]);
    
    
    subplot(2,2,4); loglog(axisX,[err_rrr(1:i);zeros(max_ind-i,1)],'b');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('Recovery error')
    xlim([512, max_ind*512]);
    
    shg;
    Frame(i+1) = getframe(fig);
    writeVideo(v,Frame(i+1));
end


%% movie

close(v);
[h, w, p] = size(Frame(1).cdata);  % use 1st frame to get dimensions
hf = figure;
%resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
set(hf, 'position', [150 50 600 600]);
axis off
fps = 10;
%movie(Frame,1,fps);
mplay(Frame)

