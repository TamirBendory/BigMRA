clear; close all; clc;

%% Defining the problem

L = 21; % length of signal
k = 10;  %2*500000; % # of signal's repetitions (maximal number)
sigma = 0.1;  % noise level2
N = 500;

%% Generating data
tic
x = randn(L, 1); 
[y, yc, ind] = gen_data(x,N,k,sigma,2*L);
snr = norm(yc)^2/norm(y-yc)^2; % The problem's SNR
k_eff = length(ind); % The actual nunber of signal's repetitions
fprintf('Measurement length  = %d \n',N);
fprintf('The measurement contains %d repetitions of the underlying signal (from %d)\n',k_eff,k);
fprintf('SNR = %.4f \n',snr);
fprintf('Generating data time = %.2f [sec] \n',toc);

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
err = zeros(num_iter+1,1);

% Initialize
% xest(:, 1) = x  + 0.1*randn(L, 1);
% xest(:, 1) = randn(L, 1);
xest(:, 1) = zeros(L, 1);

err(1) = norm(x - xest(:,1))/norm(x);
w = zeros(N, 1); 
sigma_em = sigma;

for iter = 1:num_iter
    X = repmat(xest(:, iter), 1, N)';
    S = (-0.5/(sigma_em^2))*sum((X-y_mat).^2, 2);
    S = S - max(S(:));
    w = exp(S);
    w = w/sum(w);
    xest(:,iter+1) = w'*y_mat;
    err(iter+1) = norm(x - xest(:,iter+1))/norm(x);
    if norm(xest(:, iter+1) - xest(:, iter)) < 1e-15 * norm(xest(:, iter)) % relative difference
        break
    end
end

%%
fprintf('# of iterations %u (out of %u)\n',iter,num_iter);
figure;
subplot(1,2,1);
plot(err(1:iter+1)); title('error');
subplot(1,2,2);
plot([x, xest(:, iter+1)]); legend('truth', 'estimated');