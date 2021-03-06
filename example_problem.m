clear; close all; clc;

%% Defining the problem

L = 21; % length of signal
k = 2; %5e6; %2*500000; % # of signal's repetitions (maximal number)
window_size = 4*L; 
Nfactor = 2; % Sparsity factor, should be ~6
N = window_size*k*Nfactor; % # of measurements
overlapping_factor = 1; % windows are overlapped by window_size/overlapping_factor
%% Generating data
tic
%x = randn(L,1); 
x = ones(L,1);
[y,yc, ind] = gen_data(x,N,k,0,window_size);

sigma_vec = [0.1 , 0.7, 3];  % noise level2

figure; 

for ii = 1:length(sigma_vec)
sigma = sigma_vec(ii);
y = yc + sigma*randn(N,1);
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

%% plotting

ln = 2;

subplot( length(sigma_vec),2,(ii-1)*2+1); title(strcat('sigma = ', num2str(sigma)));
hold on; plot(1:N,yc,'linewidth',ln); plot(1:N,y);
legend('clean data','data');
%title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k_eff),', SNR=',num2str(snr)));
axis tight
w = flipud(xcorr(x,y));
subplot(length(sigma_vec),2,(ii-1)*2+2);  plot(1:N,w(N:end),'linewidth',ln); 
if ii==1 title('correlation between x and y'); end
axis tight

end

pdf_print_code(gcf, 'example.pdf', 11)