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
PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
PSX = abs(fft2(X_zp)).^2;
%figure(8786); imagesc([PS - PSX]);
err_PS = norm(PSX(:) - PS(:))/norm(PSX(:));
fprintf('error PS = %.6g\n',err_PS);
th = 0.06;
X_N = double(rgb2gray(imread('newton.jpg')));
X_N = X_N/max(X_N(:));
X_N = imresize(X_N, [L, L]);
X_N = X_N - mean(X_N(:));
X_init = zeros(W);
X_init(1:L,1:L) = X_N;
%X_init = zeros(W); X_init(1:L,1:L) = -1 + 2*rand(L);
%X_init = X_zp + 0.1*randn(W);
[Xest_rrr, discrepancy_norm,err,err1, err2] = RRR_with_err(sqrt(PS),L,th,X_init,X);
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
if 1
lambda = 100;
niter = 10;
Xest_TV = TVL1denoise(Xest_rrr, lambda, niter);
%options.lambda = 1;
%[Xest_TV,err,tv,lalist] = perform_tv_denoising(Xest_rrr,options);
err_TV = norm(Xest_TV(:) - X(:))/norm(X(:));
fprintf('error TV = %.4g\n',err_TV);
end
%% display
figure(1); colormap gray;

%subplot(221);imagesc(X_init(1:L,1:L)); axis equal tight off; title('initial guess');
subplot(121);imagesc(X); axis equal tight off; title('target image');
subplot(122);imagesc(Xest_rrr); axis equal tight off; title('estimated RRR image');
%subplot(313);imagesc(Xest_TV); axis equal tight off; title('estimated RRR + TV image');


%% saving


%figure(56);imagesc(Xest_rrr); axis image off; colormap gray
%pdf_print_code(gcf, 'Einstein_1e4.pdf', 11)
