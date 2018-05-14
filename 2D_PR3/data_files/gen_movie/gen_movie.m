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
max_latte = 92;
max_polar = 22;
max_polar_v2 = 24;
max_chai = 42;
total_ind = max_latte+max_polar+max_chai+max_polar_v2;
axisX = 100*(1:total_ind);
err_PS = zeros(total_ind,1);
err_rrr = zeros(total_ind,1);

%% latte

for i = 1:max_latte
    % loading data
    str = strcat('../data_exp_',num2str(100*i),'_latte');
    data = load(str);
    M2 = data.M2;
    N = data.N;
    m_eff = data.m_eff;
    m_eff = m_eff(m_eff>0);
    m_total = sum(m_eff);
    n_micrographs = length(m_eff);
    if i == 1
        list2 = data.list2;
        sigma = data.sigma;
        W = data.W;
        X = data.X;
        L = data.L;
        X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
        PSX = abs(fft2(X_zp)).^2;
        % initial guess
        X_N = double(rgb2gray(imread('newton.jpg')));
        X_N = X_N/max(X_N(:));
        X_N = imresize(X_N, [L, L]);
        X_N = X_N - mean(X_N(:));
        X_init = zeros(W);
        X_init(1:L,1:L) = X_N;
        subplot(2,2,[1,2]); imagesc(X_N);
        colormap gray; axis equal tight off square
        shg;
        Frame(i) = getframe(fig);
        writeVideo(v,Frame(i));
    end
    % computing PS and plotting
    PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
    err_PS(i) = norm(PSX(:) - PS(:))/norm(PSX(:));
    %pause(.1);
    subplot(2,2,3);
    loglog(axisX,err_PS,'k');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('PS error')
    xlim([100, total_ind*100]);
    % RRR
    [X_out, discrepancy_norm,err,err1,err2] = RRR_with_err(sqrt(PS),L,[],X_init,X);
    [err_rrr(i),min_iter] = min(err);
    subplot(2,2,4); loglog(axisX,err_rrr,'b');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('Recovery error')
    xlim([100, total_ind*100]);
    subplot(2,2,[1,2]); imagesc(X_out); colormap gray; axis equal tight off square
    shg;
    if i>1
        Frame(i) = getframe(fig);
        writeVideo(v,Frame(i));
    end
end

%% polar

PS_latte = PS;
fprintf('Switching to polar\n');

for i = 1:max_polar
    % loading data
    str = strcat('../data_exp_',num2str(100*i),'_polar');
    data = load(str);
    M2 = data.M2;
    N = data.N;
    m_eff = data.m_eff;
    m_eff = m_eff(m_eff>0);
    m_total = sum(m_eff);
    n_micrographs = length(m_eff);
    % computing PS and plotting
    PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
    PS = (max_latte*PS_latte + i*PS)/(i+max_latte);
    err_PS(i+max_latte) = norm(PSX(:) - PS(:))/norm(PSX(:));
    %pause(.1);
    subplot(2,2,3);
    loglog(axisX,err_PS,'k');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('PS error')
    xlim([100, total_ind*100]);
    % RRR
    [X_out, discrepancy_norm,err,err1,err2] = RRR_with_err(sqrt(PS),L,[],X_init,X);
    [err_rrr(i+max_latte),min_iter] = min(err);
    subplot(2,2,4); loglog(axisX,err_rrr,'b');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('Recovery error')
    xlim([100, total_ind*100]);
    subplot(2,2,[1,2]); imagesc(X_out); colormap gray; axis equal tight off square
    shg;
    Frame(i+max_latte) = getframe(fig);
    writeVideo(v,Frame(i+max_latte));
    
end

%% chai

PS_polar = PS;
fprintf('Switching to chai\n');
max_ind = max_latte+max_polar;

for i = 1:max_chai
    % loading data
    str = strcat('../data_exp_',num2str(100*i),'_chai');
    data = load(str);
    M2 = data.M2;
    N = data.N;
    m_eff = data.m_eff;
    m_eff = m_eff(m_eff>0);
    m_total = sum(m_eff);
    n_micrographs = length(m_eff);
    % computing PS and plotting
    PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
    PS = (max_ind*PS_polar + i*PS)/(i+max_ind);
    err_PS(i+max_ind) = norm(PSX(:) - PS(:))/norm(PSX(:));
    %pause(.1);
    subplot(2,2,3);
    loglog(axisX,err_PS,'k');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('PS error')
    xlim([100, total_ind*100]);
    % RRR
    [X_out, discrepancy_norm,err,err1,err2] = RRR_with_err(sqrt(PS),L,[],X_init,X);
    [err_rrr(i+max_ind),min_iter] = min(err);
    subplot(2,2,4); loglog(axisX,err_rrr,'b');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('Recovery error')
    xlim([100, total_ind*100]);
    subplot(2,2,[1,2]); imagesc(X_out); colormap gray; axis equal tight off square
    shg;
    Frame(i+max_ind) = getframe(fig);
    writeVideo(v,Frame(i+max_ind));    
end

%% polar v2

PS_chai = PS;
fprintf('Switching to polar v2\n');
max_ind = max_latte+max_polar+max_chai;

for i = 1:max_polar_v2
    % loading data
    str = strcat('../data_exp_',num2str(100*i),'_polar_v2');
    data = load(str);
    M2 = data.M2;
    N = data.N;
    m_eff = data.m_eff;
    m_eff = m_eff(m_eff>0);
    m_total = sum(m_eff);
    n_micrographs = length(m_eff);
    % computing PS and plotting
    PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
    PS = (max_ind*PS_chai + i*PS)/(i+max_ind);
    err_PS(i+max_ind) = norm(PSX(:) - PS(:))/norm(PSX(:));
    %pause(.1);
    subplot(2,2,3);
    loglog(axisX,err_PS,'k');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('PS error')
    xlim([100, total_ind*100]);
    % RRR
    [X_out, discrepancy_norm,err,err1,err2] = RRR_with_err(sqrt(PS),L,[],X_init,X);
    [err_rrr(i+max_ind),min_iter] = min(err);
    subplot(2,2,4); loglog(axisX,err_rrr,'b');
    axis equal tight square
    xlabel('# micrographs')
    ylabel('Recovery error')
    xlim([100, total_ind*100]);
    subplot(2,2,[1,2]); imagesc(X_out); colormap gray; axis equal tight off square
    shg;
    Frame(i+max_ind) = getframe(fig);
    writeVideo(v,Frame(i+max_ind));    
end

%% movie

close(v);
[h, w, p] = size(Frame(1).cdata);  % use 1st frame to get dimensions
hf = figure;
%resize figure based on frame's w x h, and place at (150, 150)
%set(hf, 'position', [150 150 w h]);
set(hf, 'position', [150 50 600 600]);
axis off
fps = 4;
%movie(Frame,1,fps);
mplay(Frame)


% figure; loglog(axisX,err_rrr,'b');
% axis equal tight square
% xlabel('# micrographs')
% ylabel('Recovery error')
% xlim([100, total_ind*100]);
% pdf_print_code(gcf, 'err_rrr_per_micrograph', 11)


% figure; loglog(axisX,err_PS,'b');
% axis equal tight square
% xlabel('# micrographs')
% ylabel('PS estimation error')
% xlim([100, total_ind*100]);
% pdf_print_code(gcf, 'err_PS_per_micrograph', 11)
