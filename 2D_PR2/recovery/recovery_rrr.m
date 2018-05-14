% Script to run the RRR

merging_data_files;

fprintf('loading data\n');
fprintf('sigma = %.2g\n',sigma);
m_total = sum(m_eff);
fprintf('m total = %d\n',m_total);
X_zp = [X zeros(L, W-L) ; zeros(W-L, W)];
m_eff = m_eff(m_eff>0);
n_micrographs = length(m_eff);
fprintf('n_micrographs = %d\n',n_micrographs);
% 
% if isempty(gcp('nocreate'))
%    parpool(2, 'IdleTimeout', 240);
% end

%% RRR
% Make sure that N is the right number...
PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
PSX = abs(fft2(X_zp)).^2;
Freq_mask = ones(W);
lower_ind = 101; upper_ind = 100;
Freq_mask(lower_ind:upper_ind,lower_ind:upper_ind) = 0;
%Freq_mask = abs(fft2(X_zp))/max(max(fft2(X_zp)));
PS = PS.*Freq_mask;
%PSX_LP = PSX.*Freq_mask;
%PS = PSX; 
%figure(8786); imagesc([PS - PSX]);
err_PS = norm(PSX(:) - PS(:))/norm(PSX(:));
fprintf('error PS = %.6g\n',err_PS);

th = 0.03;
%X_N = double(rgb2gray(imread('newton.jpg')));
%X_N = X_N/max(X_N(:));
%X_N = imresize(X_N, [L, L]);
%X_init = zeros(W);
%X_init(1:L,1:L) = X_N;
X_init = zeros(W); X_init(1:L,1:L) = rand(L);
%X_init = X_zp;
%hgram = hist(X_zp(:),50);
%Xf = fft2(X_zp);
%Xf_phase = angle(Xf(:));
%hgram_phases = hist(Xf_phase,50);

X_LP = real(ifft2(fft2(X_zp).*Freq_mask));
X_LP = X_LP(1:L,1:L);
%X_LP = X;
[Xest_rrr, discrepancy_norm,err] = RRR_with_err(sqrt(PS),L,th,X_init,X_LP);

figure(43);
subplot(211); semilogy(1:length(discrepancy_norm),discrepancy_norm); axis tight;
subplot(212); semilogy(1:length(err),err); axis tight;

[err_rrr,min_iter] = min(err);

% Xest_rrr(1:L,1:L);
% err_rrr1 = norm(Xest_rrr - X,'fro')/norm(X(:));
% err_rrr2 = norm(rot90(Xest_rrr,2) - X,'fro')/norm(X(:));
% if err_rrr1<err_rrr2
%     err_rrr = err_rrr1;
% else
%     err_rrr = err_rrr2;
%     Xest_rrr = rot90(Xest_rrr,2);
% end
fprintf('error RRR = %.4g\n',min(err));

%% TV refinement
lambda = 1;
niter = 100;
Xest_TV = TVL1denoise(Xest_rrr, lambda, niter);
%options.lambda = 1;
%[Xest_TV,err,tv,lalist] = perform_tv_denoising(Xest_rrr,options);
err_TV = norm(Xest_TV(:) - X(:))/norm(X(:));
fprintf('error TV = %.4g\n',err_TV);

Xest_TV_sharpen = imsharpen(Xest_TV,'Radius',1,'Amount',2);
err_TV_sharp = norm(Xest_TV_sharpen(:) - X(:))/norm(X(:));
fprintf('error TV sharp = %.4g\n',err_TV_sharp);

%% display
figure(1); colormap gray;

%subplot(221);imagesc(X_init(1:L,1:L)); axis equal tight off; title('initial guess');
subplot(221);imagesc(X); axis equal tight off; title('target image');
subplot(222);imagesc(Xest_rrr); axis equal tight off; title('estimated RRR image');
subplot(223);imagesc(Xest_TV); axis equal tight off; title('estimated TV image');
subplot(224);imagesc(Xest_TV_sharpen); axis equal tight off; title('estimated sharepned image');

%subplot(133);imagesc(Xest_LS2);axis equal tight off;


%% LS 

N_eff = N*sqrt(n_micrographs);
if 0
X0 = zeros(W);
X0(1:L,1:L) = Xest_rrr;
opt_window = W;
max_iter = 20;
[Xest_LS1, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'conjugate_gradient',opt_window,max_iter);

Xest_LS1 = Xest_LS1(1:L,1:L);
err_LS1 = norm(Xest_LS1(:) - X(:))/norm(X(:));
fprintf('error LS1 = %.4g\n',err_LS1);
end

if 0
opt_window = L;
max_iter = 100;
%X0 = Xest_rrr;
X0  = Xest_LS1;
[Xest_LS2, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'conjugate_gradient',opt_window,max_iter);

err_LS2_2 = norm(Xest_LS2(:) - X(:))/norm(X(:));
fprintf('error LS2 = %.4g\n',err_LS2_2);

%save('RRR_output');

end


%% saving figures 
ind = 1000:1200;

if 0
    
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

end
