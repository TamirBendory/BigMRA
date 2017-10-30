clear; close all; clc;

L = 6; % length of signal
x = randn(L,1);

N = 200; F = dftmtx(N);
k = 20;
sigma = 0.1;

[y, yc, ind] = gen_data(x,N,k,L,sigma);
s = zeros(N,1); s(ind) = 1;
ys = [y; y(1:L-1)];
%% Non-local means

num_iter =4; err = zeros(num_iter+1,1);
err(1) = norm(y - yc)/norm(yc);
w = zeros(N);
yest = ys;

for iter = 1:num_iter
    yn = zeros(N+L-1,1);
    for i = 1:N
        yi = yest(i:i+L-1);
        for j = 1:N
            yj = yest(j:j+L-1);
            w(i,j) = exp(-0.5/(sigma^2)*norm(yi - yj)^2);
        end
    end
    w = w/sum(w(:))*N;
    for i = 1:N
        for j= 1:N
            yn(i:i+L-1) = yn(i:i+L-1) + w(i,j)*yest(j:j+L-1);
        end
    end
    
    yest = yn;
    err(iter+1) = norm(yc - yest(1:N))/norm(yc);
end
