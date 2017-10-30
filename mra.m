clear; close all; clc;

L = 4; % length of signal
x = randn(L,1); x = x - mean(x);

N = 10000;

k = 300;

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


%%  

x0 = randn(L,1);
x0 = x0/norm(x0)*norm(x);
xest = x0;
num_iter = 20;
err = zeros(num_iter,1);
%C = err; 
max_corr = zeros(N,1);

for iter = 1:num_iter
    xest_iter = xest;
    for i = 0:N-1
        
        [y_aligned, max_corr(i+1), ind] = align_to_reference(y( mod((i: i+L-1),N)+1), xest);
        xest_iter = xest_iter +  y_aligned*max_corr(i+1);
        
    end
    
    C(iter) = sum(max_corr);
    xest = xest_iter;
    xest = xest/norm(xest)*norm(x);
    xest = xest - mean(xest);
    
    [xest_aligned,~, ~] = align_to_reference(xest, x);
    err(iter) = norm(xest_aligned - x)/norm(x);
    
end

figure; plot(err); 
figure; hold on; plot(x); plot(xest_aligned); xlim([1,L])
figure; plot(C);  title('C');