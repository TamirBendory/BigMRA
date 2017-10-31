clear; close all; clc;

% Comment: Need to find a norm estimator!

%% Defining the problem

L = 21; % length of signal
k = 2e4; %2*500000; % # of signal's repetitions (maximal number)
sigma = 1;  % noise level2
window_size = 4*L; 
Nfactor = 6; % Sparsity factor, should be ~6
N = window_size*k*Nfactor; % # of measurements
overlapping_factor = 2; % windows are overlapped by window_size/overlapping_factor
%% Generating data
tic
x = randn(L,1); 
[y,yc, ind] = gen_data(x,N,k,sigma,window_size);
snr = norm(yc)^2/norm(y)^2; % The problem's SNR
k_eff = length(ind); % The actual nunber of signal's repetitions
fprintf('Measurement length  = %e \n',N);
fprintf('The measurement contains %.1e repetitions of the underlying signal (from %.1e)\n',k_eff,k);
fprintf('SNR = %.4f \n',snr);
fprintf('Generating data time = %.2f [sec] \n',toc);
%estimating the signal's norm from the data
normX = sqrt((norm(y)^2 - sigma^2*N)/k_eff);
Err_normX = abs(norm(x) - normX)/norm(x);
fprintf('Error of norm''s estimation = %.4f  \n',Err_normX);

%sanity check: note that yc == cconv(x,s,N);
%s = zeros(N,1); s(ind) = 1;
%norm(yc - cconv(x,s,N))

%% rearraging the matrix data

tic
y_mat = gen_data_mtx(y,window_size,overlapping_factor);
assert(size(y_mat,1)==window_size && size(y_mat,2)==N/window_size*overlapping_factor, 'Something is wrong with the data dimensions');
fprintf('Constructing data matrix time = %.2f [sec] \n',toc);     

%% invariants

if isempty(gcp('nocreate'))
            parpool(2, 'IdleTimeout', 240);
end

tic
[mean_est, P_est, B_est] = invariants_from_data(y_mat, sigma);
[z, problem] = phases_from_bispectrum_real(B_est, sign(mean_est), randn(window_size,1));
x_est = real(ifft(sqrt(P_est).*z));
xref = zeros(window_size,1); xref(1:L) = x;
[x_aligned,~, ~] = align_to_reference(x_est, xref);
x_aligned = x_aligned(1:L); 
% using norm(x) for for a while
%x_aligned = x_aligned/norm(x_aligned)*norm(x);
% correcting with norm estimation
x_aligned = x_aligned/norm(x_aligned)*normX;
fprintf('Algorithm time = %.2f [sec] \n',toc);

%% plotting

err = norm(x_aligned - x)/norm(x);
inds = ind(10) - 300; indf = inds + 600;

figure; 
subplot(211); hold on; stem(1:L,x); stem(1:L,x_aligned,'xr'); 
title(strcat('Error = ',num2str(err)));
legend('signal','estimation');
subplot(212); hold on; plot(inds:indf,yc(inds:indf),'linewidth',2); plot(inds:indf,y(inds:indf)); legend('clean data','data');
title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k_eff),', SNR=',num2str(snr)));
axis tight
    