% Script to run the RRR

load('data_exp_gpu_first_round.mat');
%load('data_exp_gpu.mat');

%merging_data_files;

fprintf('loading data\n');
fprintf('sigma = %.2g\n',sigma);
m_total = sum(m_eff);
fprintf('m total = %d\n',m_total);
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
m_eff = m_eff(m_eff>0);
n_micrographs = length(m_eff);
fprintf('n_micrographs = %d\n',n_micrographs);

if isempty(gcp('nocreate'))
    parpool(2, 'IdleTimeout', 240);
end

%% RRR
% Make sure that N is the right number...
PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
PSX = abs(fft2(X_zp)).^2;
%figure(8786); imagesc([PS - PSX]);
err_PS = norm(PSX(:) - PS(:))/norm(PSX(:));
fprintf('error PS = %.6g\n',err_PS);

dc = M1/m_total; % sum(X(:))
th = .023;

[Xest_rrr, discrepancy_norm] = RRR(sqrt(PS),L,th);
figure(43); semilogy(1:length(discrepancy_norm),discrepancy_norm); axis tight;
%figure(44); plot(1:length(discrepancy_norm),discrepancy_norm); axis tight;

Xest_rrr = Xest_rrr(1:L,1:L);
err_rrr1 = norm(Xest_rrr - X,'fro')/norm(X(:));
err_rrr2 = norm(rot90(Xest_rrr,2) - X,'fro')/norm(X(:));
if err_rrr1<err_rrr2
    err_rrr = err_rrr1;
else
    err_rrr = err_rrr2;
    Xest_rrr = rot90(Xest_rrr,2);
end
fprintf('error RRR = %.4g\n',err_rrr);


%% LS 

N_eff = N*sqrt(n_micrographs);
if 1
X0 = zeros(W);
X0(1:L,1:L) = Xest_rrr;
opt_window = W;
max_iter = 200;
[Xest_LS1, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'trust_regions',opt_window,max_iter);

Xest_LS1 = Xest_LS1(1:L,1:L);
err_LS1 = norm(Xest_LS1(:) - X(:))/norm(X(:));
fprintf('error LS1 = %.4g\n',err_LS1);
end

opt_window = L;
max_iter = 1000;
%X0 = Xest_rrr;
X0  = Xest_LS1;
[Xest_LS2, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'trust_regions',opt_window,max_iter);

err_LS2_2 = norm(Xest_LS2(:) - X(:))/norm(X(:));
fprintf('error LS2 = %.4g\n',err_LS2_2);

save('RRR_output');


%% display
figure(1); colormap gray;
subplot(131);imagesc(X); axis equal tight off;
subplot(132);imagesc(Xest_LS1); axis equal tight off;
subplot(133);imagesc(Xest_LS2);axis equal tight off;

%% saving figures 
ind = 1000:1200;

figure(56);imagesc(Y_clean(ind,ind)); axis image off; colormap gray
pdf_print_code(gcf, 'data2D_clean.pdf', 11)

sigma = 0.2;
figure(57); imagesc(sigma*randn(length(ind))+ Y_clean(ind,ind)); axis image off; colormap gray
pdf_print_code(gcf, 'data2D_noisy_02.pdf', 11)

figure(58); imagesc(Y_obs(ind,ind)); axis image off; colormap gray
pdf_print_code(gcf, 'data2D_noisy_1.pdf', 11)

figure(59); imagesc(X); axis image off; colormap gray
pdf_print_code(gcf, 'signal2D_clean.pdf', 11)

figure(60); imagesc(Xest_rrr); axis image off; colormap gray
pdf_print_code(gcf, 'signal2D_RRR.pdf', 11)

figure(61); imagesc(Xest_LS2); axis image off; colormap gray
pdf_print_code(gcf, 'signal2D_LS_intermediate.pdf', 11)


