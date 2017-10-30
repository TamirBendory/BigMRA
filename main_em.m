clear; close all; clc;

L = 10; % length of signal
x = ones(L,1); %x = x/norm(x);
%x = randn(L,1); x = x/norm(x);


N = 4000000;
k = 8000;
sigma = 4;

[y,yc, ind] = gen_data(x,N,k,L,sigma);
snr = norm(yc)^2/norm(y)^2;
k_eff = length(ind);
fprintf('The data contains %d repetations of the underlying signal\n',k_eff);
%s = zeros(N,1); s(ind) = 1;
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
num_iter = 1000;
xest = zeros(L,num_iter+1);
xest(:,1) = randn(L,1);
err = zeros(num_iter+1,1);
err(1) = norm(x - xest(:,1))/norm(x);
w = zeros(N,1);
sigma_em = 4*sigma;

for iter = 1:num_iter
    X = repmat(xest(:,iter),1,N)';
    w = exp(-0.5/(sigma_em^2)*sum((X-y_mat).^2,2));
%     for i = 1:N
%         % here I can choose different possibilities
%         w(i) = exp(-0.5/(sigma_em^2)*norm(xest(:,iter)'-y_mat(i,:))^2);
%         %w(i) = 1./(norm(xest(:,iter)-y_mat(i,:))^2 + sigma_em);       
%     end
    w = w/sum(w);
    xest(:,iter+1) = w'*y_mat; xest(:,iter+1) = xest(:,iter+1)/norm(xest(:,iter+1))*norm(x);
    err(iter+1) = norm(x - xest(:,iter+1))/norm(x);
    if abs(err(iter+1) - err(iter))<10^-6
        break
    end
end

figure; 
subplot(311); hold on; stem(x); stem(xest(:,iter),'xr'); legend('signal','estimation');
subplot(312); hold on; plot(yc(1:1000),'linewidth',4); plot(y(1:1000)); legend('clean data','data');
title(strcat('N =', num2str(N),', L=',num2str(L), ', K=',num2str(k_eff),', SNR=',num2str(snr)));
subplot(313);plot(1:iter, err(1:iter)); title(strcat('err =', num2str(err(iter+1))));

