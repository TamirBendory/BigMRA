clear; close all; clc;

%% Defining the problem

L = 21; % length of signal
k = 1e6; % # of signal's repetitions (maximal number)
sigma = 2.5;  % noise level2
W = 8*L; 
Nfactor = 6; % Sparsity factor, should be ~6
N = W*k*Nfactor; % # of measurements
overlapping_factor = 1; % windows are overlapped by window_size/overlapping_factor
%% Generating data
tic
x = randn(L,1);
%x = ones(L,1);
[y,yc, ind] = gen_data(x,N,k,sigma,W);
snr = norm(yc)^2/norm(y-yc)^2; % The problem's SNR
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
if overlapping_factor == 1
y_mat = reshape(y, W,N/W);
else
y_mat = gen_data_mtx(y,W,overlapping_factor);
fprintf('Constructing data matrix time = %.2f [sec] \n',toc);     
end
assert(size(y_mat,1)==W && size(y_mat,2)==N/W*overlapping_factor, 'Something is wrong with the data dimensions');

%% invariants

if isempty(gcp('nocreate'))
            parpool(2, 'IdleTimeout', 240);
end

tic
[mean_est, P_est, B_est] = invariants_from_data(y_mat, sigma);
[z, problem] = phases_from_bispectrum_real(B_est, sign(mean_est), randn(W,1));
x_est = real(ifft(sqrt(P_est).*z));

%% automatic alignment

x_aligned = auto_alignment(x_est,L,0,x);
x_aligned = x_aligned/norm(x_aligned)*normX;
fprintf('Algorithm time = %.2f [sec] \n',toc);

%% plotting

err = norm(x_aligned - x)/norm(x)
lag = 600;
inds = ind(10) - lag/2; indf = inds + lag;
yc_ind = yc(inds:indf);
y_ind = y(inds:indf);

save('x.mat','x');
save('x_aligned.mat','x_aligned');
save('err.mat','err');
save('y_ind.mat','y_ind');
save('yc_ind.mat','yc_ind');

%% 

figure; 
subplot(311); hold on; stem(1:L,x); stem(1:L,x_aligned,'xr'); 
title(strcat('Error = ',num2str(err)));
legend('signal','estimation');
subplot(312); hold on; 
%plot(inds:indf,yc_ind,'linewidth',2); plot(inds:indf,y_ind); 
plot(1:601,yc_ind,'linewidth',2); plot(1:601,y_ind); 
legend('clean data','data');
title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k),', SNR=',num2str(snr)));
axis tight
w = flipud(xcorr(x,y_ind));
subplot(313);  
%plot(inds:indf,w(lag+1:end)); 
plot(1:601,w(lag+1:end)); 
title('correlation between x and y');
axis tight
   