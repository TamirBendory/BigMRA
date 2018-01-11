function [y, y_clean, ind, class] = gen_data2D(x, N, m, sigma, W)
% Generate data following the Big-MRA model - the 2D version
%
% Input:
%
% x - the K signals to be hidden in noise (LxK matrix)
% N - total length of the observation
% m - vector of length K such that x(:, k) is repeated (nearly) m(k) times
% sigma - noise level
% W - length of intended sliding window: only used to ensure proper separation
%
% Output:
%
% y - noisy measurement of length N
% y_clean - clean (noiseless) measurement
% ind - locations of x in y_clean
% class - for each location in ind, specifies which signal (1..K) is there.

y_clean = zeros(N);
% currently, we work only for the homogenous case
L = size(x,1);
ind = randi(N-L+1, 1000*m, 2);
%ind = sortrows(ind);
indn = ind(1,:);

% removing adjacent signals (need to be rewritten)
sep = W;
for i = 1:length(ind)-1
    
    if size(indn,1) == m
        break;
    end
    
    %if min(min(min(mod(repmat(ind(i+1,:),size(indn,1),1) - indn,N)))...
           % ,min(min(mod((indn - repmat(ind(i+1,:),size(indn,1),1)),N))))>sep 
%     if min(min(min(mod( bsxfun(@minus, indn, ind(i+1,:)),N)))...
%             ,min(min(mod(bsxfun(@minus, -indn, ind(i+1,:)),N))))>sep         
%         indn = [indn; ind(i+1,:)];
%     end
    ind_temp =  abs(bsxfun(@minus,indn,ind(i+1,:)));
    if min(ind_temp(:))>sep         
        indn = [indn; ind(i+1,:)];
    end
end
%     if ind(i+1,1) - ind(i,1) > sep
%         indn = [indn; ind(i+1)];
%     end
%end
ind = indn;
% drawing classes
%class = randi([1,K], length(ind), 1);
class = [];
% generating noiseless measurement(depends on the window length)
for i = 1:size(ind,1)
    y_clean( ind(i,1) : ind(i,1)+L-1,ind(i,2) : ind(i,2)+L-1 ) = y_clean( ind(i,1) : ind(i,1)+L-1,ind(i,2) : ind(i,2)+L-1 )...
        +x ; %x(:, class(i));
end

% add noise
y = y_clean + sigma*randn(N);

end
