clear; close all; clc;
dbstop if error
L = 10; % length of signal
%x = randn(L,1); x = x/norm(x);
x = ones(L,1); %x = x/norm(x);

N = 50000;
k = 100;
sigma = .1;


[y ind] = gen_data(x,N,k,L,sigma);


%% running over indices
max_iter = 1000;
x_est = randn(L,1);
ys = [y ; y(1:L)];
w = zeros(N,max_iter);
err = zeros(max_iter,1);

for iter = 1:max_iter
    
    for i = 1:N
        w(i,iter) = exp(-norm(ys(i:i+L-1) - x_est)/2/(sigma^2));
    end
    
    w(:,iter) = w(:,iter)/sum(w(:,iter)); %normalization
    
    x_new = zeros(L,1);
    for i = 1:N
        x_new = x_new + ys(i:i+L-1)*w(i,iter);
    end
    x_est = x_new/norm(x_new)*norm(x); 
    %x_est = x_new/norm(x_new)*norm(x);
    err(iter) = norm(x - x_est);% /norm(x);
end

figure; plot(err)
%%

figure; hold on;
stem(x); stem(x_est,'*r');
legend('signal','recovered');