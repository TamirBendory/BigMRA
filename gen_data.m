function [y, yc, ind, class] = gen_data(x,N,m,sigma,W,K)

% generating the problem measurements
% input:
% x - the underlying signal
% N - length of measurement
% m - number of repetitions
% sigma - noise level
% W - the size of analysis window
% K - # of different signals

% output:
% y - noisy measurement
% yc - clean (noiseless) measurement
% ind - location of x in yc

yc = zeros(N,1);
L = length(x);
ind = randi([L,N-L], m, 1); ind  = sort(ind);
indn = ind(1);

% removing adjacent signals (need to be rewritten)

for i = 1:length(ind)-1
    if ind(i+1) - ind(i) > 2*W
        indn = [indn; ind(i+1)];
    end
end
ind = indn;
% drawing classes
class = randi([1,K], length(ind), 1);

% generating noiseless measurement(depends on the window length)
for i = 1:length(ind)
    yc( ind(i): ind(i)+L-1) = yc( ind(i): ind(i)+L-1) +  x(:,class(i));
end

% add noise
y = yc + sigma*randn(N,1);

end