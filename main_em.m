clear; close all; clc;

% Need to find norm estimator!

L = 11; % length of signal
B = 1; % bandwidth
x = fft(randn(L,1)); x(B+2:L-B) = 0; 
x = real(ifft(x));

%x = zeros(L,1); x(1) = 1;
%x = ones(L,1); %x = x/norm(x);
%x = randn(L,1); x = x/norm(x);

N = 4000;
k = 200;
sigma = 0.01;

[y,yc, ind] = gen_data(x,N,k,L,sigma);
snr = norm(yc)^2/norm(y)^2;
k_eff = length(ind);
fprintf('The data contains %d repetations of the underlying signal\n',k_eff);
s = zeros(N,1); s(ind) = 1;
% note that y == cconv(x,s,N);
% sanity check:
%norm(yc - cconv(x,s,N))
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
xest(:,1) = fft(randn(L,1)); xest(B+2:L-B,1) = 0; 
xest(:,1) = real(ifft(xest(:,1)));
% xest(:,1) = randn(L,1);
% xest(:,1) = xest(:,1)/norm(xest(:,1))*norm(x);
err = zeros(num_iter+1,1);
err(1) = norm(x - xest(:,1))/norm(x);
w = zeros(N,1);
sigma_em = 100*sigma;

for iter = 1:num_iter
    X = repmat(xest(:,iter),1,N)';
    w = exp(-0.5/(sigma_em^2)*sum((X-y_mat).^2,2));
    w = w/sum(w);% w = circshift(w,(L-1)/2-2);
    %Maxw = max(w);
    %w(w<0.01*Maxw) = 0;
    xest(:,iter+1) = w'*y_mat; xest(:,iter+1) = xest(:,iter+1)/norm(xest(:,iter+1))*norm(x);
    err(iter+1) = norm(x - xest(:,iter+1))/norm(x);
    if abs(err(iter+1) - err(iter))<10^-8
        break
    end
end


%% plotting

figure; 
subplot(411); hold on; stem(x); stem(xest(:,iter),'xr'); legend('signal','estimation');
subplot(412); hold on; plot(yc(1:500),'linewidth',2); plot(y(1:500)); legend('clean data','data');
title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k_eff),', SNR=',num2str(snr)));
subplot(413);plot(1:iter, err(1:iter)); title(strcat('err =', num2str(err(iter+1))));
subplot(414);plot(1:500,w(1:500)); title('Weights of last iteration');

