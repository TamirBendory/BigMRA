clear; close all; clc;

L = 4; % length of signal
x = randn(L,1); x = x - mean(x);

N = 10000;

k = 300;

sigma = 0;

y = sigma*randn(N,1);
s = randi(N, k, 1);

s  = sort(s);
indn = s(1);
for i = 1:length(s)-1
    if s(i+1) - s(i) > 2*L
        indn = [indn; s(i+1)];
    end
end

s = indn;

for i = 1:length(s)
    y( s(i): s(i)+L-1) = y( s(i): s(i)+L-1) +  x;
end


%%  

x0 = randn(L,1);
%x0 = x0/norm(x0)*norm(x);
x0 = [x0; zeros(N-L,1)];
xest = x0;
xfest = fft(x0);

num_iter = 20;
err = zeros(num_iter,1);
%C = err; 
max_corr = zeros(N,1);
yf = fft(y);
% recall that we are trying to solve y = x \ast s +n;

for iter = 1:num_iter
    sf_est = wiener(xfest,yf,sigma);
   
    xfest = wiener(sf_est,yf,sigma);
    %[xest_aligned,~, ~] = align_to_reference(xest, x);
    %err(iter) = norm(xest_aligned - x)/norm(x);
    
end

figure; plot(err); 
figure; hold on; plot(x); plot(xest_aligned); xlim([1,L])
figure; plot(C);  title('C');