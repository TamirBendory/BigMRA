clear; close all; clc;

L = 4; % length of signal
x = randn(L,1); x = x - mean(x);

N = 1000;
k = 30;

sigma = 0.01;

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

xest = randn(L,1);
num_iter = 10000;
err = zeros(num_iter,1); C = err;  C_candidate = err;
max_corr = zeros(2*L+1,1);
max_corr_candidate = max_corr;

for iter = 1:num_iter
    ind = randi(L, 1);
    xest_candidate = xest;
    xest_candidate(ind) = randn(1);
    
   run_ind = 0;
   
    for i = ind-L:ind+L
        
        run_ind = run_ind+1;
        [~, max_corr(run_ind), ~] = align_to_reference(y( mod((i: i+L-1),N)+1), xest_candidate);
        [~, max_corr_candidate(run_ind), ~] = align_to_reference(y( mod((i: i+L-1),N)+1), xest);
        
%        xest_candidate = xest_candidate +  y_aligned*max_corr(i+1);
        
    end
    
    C(iter) = sum(max_corr);
    C_candidate(iter) = sum(max_corr_candidate);
    
    w = abs(C(iter)./ (C(iter) + C_candidate(iter)));
    
    if rand(1) > w
        xest = xest_candidate;
    end
    %xest = xest_iter;
    xest = xest/norm(xest)*norm(x);
    xest = xest - mean(xest);
    
    [xest_aligned,~, ~] = align_to_reference(xest, x);
    err(iter) = norm(xest_aligned - x)/norm(x);
    
end

figure; plot(err); 
figure; hold on; plot(x); plot(xest_aligned); xlim([1,L])
figure; plot(C);  title('C');