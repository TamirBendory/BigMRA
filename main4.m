clear; close all; clc;

L = 6; % length of signal
x = randn(L,1);

N = 200; F = dftmtx(N);
k = 20;
sigma = 0;

[y ind] = gen_data(x,N,k,L,sigma);
s = zeros(N,1); s(ind) = 1;
%% rearraging the matrix data
%y_matrix = zeros(N,L);
y_stretch = [y ; y(1:L-1)];
for i = 1:2*L:N
    y_matrix(i,:) = y_stretch(i:i+L-1);
end
clear y_stretch;

%%
yf = fft(y);
max_iter = 40; errX = zeros(max_iter,1); errS = zeros(max_iter,1);
%s_est =  rand(N,1); %rand(N,1);
warning('off');
x_est = zeros(N,1);
x_est(1:L) = randn(L,1);
lambda = 10;

for i = 1:max_iter
    
 %   sf = fft(s_est);
   
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
subplot(211); plot(errX);  title ('Error X');
subplot(212); plot(errS);  title ('Error S');
