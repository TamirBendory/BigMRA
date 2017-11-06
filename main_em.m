clear; close all; clc;

%% Defining the problem

L = 21; % length of signal
k = 1e4; %5e6; %2*500000; % # of signal's repetitions (maximal number)
sigma = 0.1;  % noise level2
window_size = 4*L; 
Nfactor = 6; % Sparsity factor, should be ~6
N = window_size*k*Nfactor; % # of measurements

%% Generating data
tic
B = 2; % bandwidth
x = fft(randn(L,1)); x(B+2:L-B) = 0; 
x = real(ifft(x)); % generating low-pass signal
[y,yc, ind] = gen_data(x,N,k,sigma,window_size);
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

%% rearraging the matrix data

y_stretch = [y ; y(1:L-1)];
y_mat = zeros(N,L);
for i = 1:N
    y_mat(i,:) = y_stretch(i:i+L-1);
end
clear y_stretch;

%% EM iterations
num_iter = 100;
xest = zeros(L,num_iter+1);
%xest(:,1) = x; %intial guess
xest(:,1) = fft(randn(L,1)); xest(B+2:L-B,1) = 0; 
xest(:,1) = real(ifft(xest(:,1)));

err = zeros(num_iter+1,1); err(1) = norm(x - xest(:,1))/norm(x);
w = zeros(N,1); sigma_em = sigma;

for iter = 1:num_iter
    X = repmat(xest(:,iter),1,N)';
    S = -0.5/(sigma_em^2)*sum((X-y_mat).^2,2);
    S = S - max(S(:));
    w = exp(S);
    w = w/sum(w);
    xest(:,iter+1) = w'*y_mat; 
    xest(:,iter+1) = xest(:,iter+1)/norm(xest(:,iter+1))*norm(x); % plugging in norm(x)
    err(iter+1) = norm(x - xest(:,iter+1))/norm(x);
    if norm(xest(:, iter+1) - xest(:, iter)) < 1e-10 * norm(xest(:, iter)) % relative difference
        break
    end
end

fprintf('# of iterations %u (out of %u)\n',iter,num_iter);

%% plotting

lag = 600;
inds = ind(10) - lag/2; indf = inds + lag;

figure; 
subplot(411); %title(strcat('Error = ',num2str(err)));
hold on; stem(1:L,x); stem(1:L,xest(:, iter+1),'xr'); 
legend('signal','estimation');
axis tight;

subplot(412); hold on; plot(inds:indf,yc(inds:indf),'linewidth',2); plot(inds:indf,y(inds:indf)); legend('clean data','data');
title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k_eff),', SNR=',num2str(snr)));
axis tight

subplot(413); semilogy(inds:indf,w(inds:indf),'linewidth',2);
% set(gca, 'YScale', 'log');
title('weights (log10 scale)');
axis tight

subplot(414); hold on; plot(err(1:iter+1),'linewidth',2); 
title('Error'); xlabel('iteration'); ylabel('error')
axis tight

fprintf('max w: %g (if close to 1, suggests only one segment of the micrograph counts, which isn''t good)\n', max(w));