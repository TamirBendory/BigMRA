function [y, yc, ind] = gen_data(x,N,k,sigma,window_size) 

% generating the problem measurements
% input:
% x - the underlying signal
% N - length of measurement
% k - number of repetitions
% sigma - noise level
% window_size - the size of analysis window
% output:
% y - noisy measurement
% yc - clean (noiseless) measurement
% ind - location of x in yc

yc = zeros(N,1);
L = length(x);
ind = randi([L,N-L], k, 1); ind  = sort(ind);
indn = ind(1);

% removing adjacent signals (need to be rewritten)
for i = 1:length(ind)-1
    if ind(i+1) - ind(i) > window_size
        indn = [indn; ind(i+1)];
    end
end
ind = indn;

% generating noiseless measurement
for i = 1:length(ind)
    yc( ind(i): ind(i)+L-1) = yc( ind(i): ind(i)+L-1) +  x;
end

% add noise
y = yc + sigma*randn(N,1);

end