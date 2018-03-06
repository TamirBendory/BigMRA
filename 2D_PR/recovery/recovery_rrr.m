% Script to run the RRR

%load('data_exp_gpu_first_round.mat');
%load('data_exp_gpu.mat');

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

%% RRR
% Make sure that N is the right number...
PS = M2PS_partial(M2,list2,m_total,sigma,n_micrographs*N^2,W);
PSX = abs(fft2(X_zp)).^2;
%figure(8786); imagesc([PS - PSX]);
err_PS = norm(PSX(:) - PS(:))/norm(PSX(:));
fprintf('error PS = %.6g\n',err_PS);

dc = M1/m_total; % sum(X(:))
th = .017;
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

%opt_window = W;
%max_iter = 100;
%[Xest_LS1, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, Xest_rrr,'trust_regions',opt_window,max_iter);

%Xest_LS1 = Xest_LS1(1:L,1:L);
%err_LS1 = norm(Xest_LS1(:) - X(:))/norm(X(:));
%fprintf('error LS1 = %.4g\n',err_LS1);

opt_window = L;
max_iter = 500;
X0 = Xest_rrr;
[Xest_LS2, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'trust_regions',opt_window,max_iter);

err_LS2_2 = norm(Xest_LS2(:) - X(:))/norm(X(:));
fprintf('error LS2 = %.4g\n',err_LS2_2);

save('RRR_output');


%% display
figure(1); colormap gray;
subplot(131);imagesc(X); axis equal tight off;
subplot(132);imagesc(Xest_rrr); axis equal tight off;
subplot(133);imagesc(Xest_LS2);axis equal tight off;



%% 

% 
% Mask = zeros(W); 
% Mask(1:L,1:L) = 1;
% 
% 
% 
% X0 = Xest_rrr;
% 
% num_epocs = 2;
% err_epoc = zeros(num_epocs,2);
% obj_epoc = zeros(num_epocs,1);
% 
% for i = 1:num_epocs
% 
% [X_est, problem,stats] = least_squares_2D(M1, M2, M3, W, sigma, N_eff, L, m_total, list2, list3, X0,'trust_regions');
% X_est = align_to_reference_2D(X_est, X_zp);
% fprintf('epoc = %d\n',i);
% err_epoc(i,1) = norm(X_zp(:) -  X_est(:))/norm(X(:));
% fprintf('err before proj = %.8g\n',err_epoc(i,1));
% obj_epoc(i) = stats.cost;
% X0 = X_est;
% X0 = X0.*Mask;
% err_epoc(i,2) = norm(X_zp(:) -  X0(:))/norm(X(:));
% fprintf('err after proj = %.8g\n',err_epoc(i,2));
% fprintf('saving\n');
% save('LS_output');
% end
% 
