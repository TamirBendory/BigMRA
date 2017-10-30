clear; close all; clc;

L = 4; % length of signal
%x = randn(L,1); 
x = randn(L,1);
x = x - mean(x);
N = 1000;
x_comp = zeros(N,1); x_comp(1:L) = x;

k = 30;

sigma = 0;

y = sigma*randn(N,1);
ind = randi(N, k, 1);

ind  = sort(ind);
indn = ind(1);
for i = 1:length(ind)-1
    if ind(i+1) - ind(i) > 2*L
        indn = [indn; ind(i+1)];
    end
end

ind = indn;

for i = 1:length(ind)
    y( ind(i): ind(i)+L-1) = y( ind(i): ind(i)+L-1) +  x;
end

%% Algorithm 

% We are trying to solve Ds = y, where s\in{0,1}. Since D is circulant, we
% can write F^*SFs = y, or SFs = Fy, where S is a diagonal matrix. We then
% aim to solve the LS problem min_{S,s} ||SFs - Fy||_2^2. Let Svec be the
% vectorization of S. Then, given s, Svec = hy/hs, where hy stands for the Fourier transform of y. 
% Given S, the hs = hy/Svec. But here we will need to take into account the
% fact that s\in{0,1}.

hy = fft(y);
num_iter = 20;
err_s = zeros(num_iter,1); err_x = err_s;
true_sup = zeros(N,1); true_sup(ind) = 1;
%s = rand(N,1); 
s = true_sup; 
eps = 10^(-10);

for i = 1:num_iter
    
    hs = fft(s);
    Svec = my_pinv(hy,hs,eps);
    hs = my_pinv(hy,Svec,eps);
    s = sqrt(N)*real(ifft(s)); 
    s(s>0.5) = 1;
    s(s<0.5) = 0;
    
    err_s(i) = norm(s - true_sup)/norm(true_sup);
    [xest,~, ~] = align_to_reference(real(ifft(Svec)), x_comp);
    err_x(i) = norm(x_comp - xest)/norm(x_comp);
    
end

%% 
% 
% C = err; 
% max_corr = zeros(N,1);
% 
% for iter = 1:num_iter
%     xest_iter = xest;
%     for i = 0:N-1
%         
%         [y_aligned, max_corr(i+1), ind] = align_to_reference(y( mod((i: i+L-1),N)+1), xest);
%         xest_iter = xest_iter +  y_aligned*max_corr(i+1);
%         
%     end
%     
%     C(iter) = sum(max_corr);
%     xest = xest_iter;
%     xest = xest/norm(xest)*norm(x);
%     xest = xest - mean(xest);
%     
%     [xest_aligned,~, ~] = align_to_reference(xest, x);
%     err(iter) = norm(xest_aligned - x)/norm(x);
%     
% end
% 
% figure; plot(err); 
% figure; hold on; plot(x); plot(xest_aligned); xlim([1,L])
% figure; plot(C);  title('C');