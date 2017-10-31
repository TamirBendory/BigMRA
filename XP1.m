clear; close all;
clc;


L = 7; % length of signal
window_size = 6*L;
Nfactor = 6;
N = window_size*k*Nfactor;%4000;
k = 100000;

sigma_vec = logspace(-3,0.5,8);
num_iter = 20;
snr = zeros(length(sigma_vec), num_iter);
err = snr;

for s = 1:length(sigma_vec)
    sigma = sigma_vec(s)
    for iter = 1:num_iter
        
        x = randn(L,1); %x = x/norm(x);
        
        %sigma = 2;
        
        [y,yc, ind] = gen_data(x,N,k,L,sigma,window_size);
        snr(s,iter) = norm(yc)^2/norm(y)^2;
        k_eff = length(ind);
        %fprintf('The data contains %d repetations of the underlying signal\n',k_eff);
        %s = zeros(N,1); s(ind) = 1;
        % note that y == cconv(x,s,N);
        % sanity check:
        %norm(yc - cconv(x,s,N))
        %% rearraging the matrix data
        
        y_stretch = [y ; y(1:window_size-1)];
        Nw = N/window_size;
        y_mat = zeros(window_size,Nw);
        for i = 0:Nw-1
            y_mat(:,i+1) = y_stretch(i*window_size+1:(i+1)*window_size);
        end
        clear y_stretch;
        
        if isempty(gcp('nocreate'))
            parpool(2, 'IdleTimeout', 240);
        end
        
        
        %% invariants
        
        [mean_est, P_est, B_est] = invariants_from_data(y_mat, sigma);
        [z, problem] = phases_from_bispectrum_real(B_est, sign(mean_est), randn(window_size,1));
        
        x_est = real(ifft(sqrt(P_est).*z));
        xref = zeros(window_size,1); xref(1:L) = x;
        
        [x_aligned,~, ~] = align_to_reference(x_est, xref);
        
        x_aligned = x_aligned(1:L); x_aligned = x_aligned/norm(x_aligned)*norm(x);
        
        %% plotting
        
        err(s,iter) = norm(x_aligned - x)/norm(x);
        
        save('err','err')
    end
end
%
% inds = ind(4) - 200;
% indf = inds + 400;
%
%
% figure;
% subplot(211); hold on; stem(x); stem(x_aligned,'xr'); title(strcat('Error = ',num2str(err)));
% legend('signal','estimation');
% subplot(212); hold on; plot(inds:indf,yc(inds:indf),'linewidth',2); plot(inds:indf,y(inds:indf)); legend('clean data','data');
% title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k_eff),', SNR=',num2str(snr)));

