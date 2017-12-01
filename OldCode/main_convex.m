clear; close all; clc;

L = 6; % length of signal
x = randn(L,1);

N = 200; F = dftmtx(N);
k = 20;
sigma = 1;

[y,yc,  ind] = gen_data(x,N,k,L,sigma);
s = zeros(N,1); s(ind) = 1;

%%
yf = fft(y);
max_iter = 20; errX = zeros(max_iter,1); errS = zeros(max_iter,1);
warning('off');
x_est = zeros(N,1);
x_est(1:L) = x + randn(L,1); %randn(L,1);
lambda = 2;

for i = 1:max_iter
    
    xf = fft(x_est);
    
        %estimating s given x
        cvx_begin
        variable s_est(N)
        minimize norm( yf - (F*s_est).*xf ,2) + lambda*norm(s_est,1)
        %subject to
        %yf  == (F*s_est).*xf;
        s_est <= 1;
        s_est >= 0;
        cvx_end
%         tau = 0.1;
%         s_est(s_est>0.5) = s_est(s_est>0.5) + tau;
%         s_est(s_est<0.5) = s_est(s_est<0.5) - tau;
%     
    errS(i) = norm(s_est - s)/norm(s);
  
    sf = fft(s_est);
    
     %estimating x given s
    cvx_begin
    variable x_est(N)
    minimize norm( yf - (F*x_est).*sf ,2)
    subject to
    x_est(L+1:end) == 0;
    cvx_end
    
    errX(i) = norm(x - x_est(1:L))/norm(x);
    
    
end

figure;
subplot(311); hold on; stem(x); stem(x_est(1:L)); 
subplot(312); plot(errX);  title ('Error X');
subplot(313); plot(errS);  title ('Error S');
