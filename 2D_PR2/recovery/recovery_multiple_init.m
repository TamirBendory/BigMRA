% Script to run the LS from multiple initializations

%load('data_exp_gpu_1500.mat');

merging_data_files;

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

%% PS estimations

PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
PSX = abs(fft2(X_zp)).^2;
%figure(8786); imagesc([PS - PSX]);
err_PS = norm(PSX(:) - PS(:))/norm(PSX(:));
fprintf('error PS = %.6g\n',err_PS);

%% LS
N_eff = N*sqrt(n_micrographs);
max_iter_1 = 20;
max_iter_2 = 100;
num_init = 100;
err_LS = zeros(num_init,1);

for i = 1:num_init
    X0  = rand(L);
    opt_window = W;
    [Xest_LS_W, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'conjugate_gradient',opt_window,max_iter_1) ; 
    opt_window = L;
    [Xest_LS, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, Xest_LS_W(1:L,1:L),'conjugate_gradient',opt_window,max_iter_2);
    err1 = norm(Xest_LS - X,'fro')/norm(X(:));
    err2 = norm(rot90(Xest_LS,2) - X,'fro')/norm(X(:));
    if err1<err2
        err_LS(i) = err1;
    else
        err_LS(i) = err2;
        Xest_LS = rot90(Xest_LS,2);
    end
    fprintf('iter = %g, error = %.4g\n',i,err_LS(i));
    if i==1 
                Xest = Xest_LS;
    elseif  err_LS(i)<min(err_LS(1:i-1))
        Xest = Xest_LS;
    end
end
    %% display
    figure(1); colormap gray;
    subplot(121);imagesc(X); axis equal tight off;
    subplot(122);imagesc(Xest); axis equal tight off;
    
    
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
