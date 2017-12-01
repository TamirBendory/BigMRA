clear; %clc; close all;

L = 10; 
x = randn(L,1);
hx = fft(x); hx = hx./abs(hx); x = real(ifft(hx));

N = 5000000;
y = zeros(L,N);
h = zeros(L,N);
sigma_vec = rand(L,1);
%sigma_vec = linspace(1/L,1,L)';
Y = zeros(L);

for n = 1:N
    
    h(:,n) = sigma_vec.*randn(L,1);
    %y(:,n) = (cconv(h(:,n),circshift(flipud(x),1),L));
    y(:,n) = cconv(h(:,n),x,L);
    %y(:,n) = circshift(flipud(y(:,n)),1);
    Y = Y + y(:,n)*y(:,n)';

end

Y = Y/N;
[U V] = eig(Y);
x_est = U(:,3);
hx_est = fft(x_est); hx_est = hx_est./abs(hx_est);
x_est = real(ifft(hx_est));
x_est = x_est*sign(mean(x))*sign(mean(x_est));

%x_est = flipud(x_est);
x_est = align_to_reference(x_est,x);

err = norm(x - x_est)/norm(x);
disp(err);

%figure; hold on;
%plot(x); plot(x_est);