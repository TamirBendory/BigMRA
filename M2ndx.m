function M = M2ndx(x,W,sigma,m,N)

% x - signal
% W - window's length
% sigma - noise level
% m - # signals' repetitions
% N - measurement length
% M1 - the first moment (I am considering the real one; not the estimated)

L = length(x); % signal's length
K = length(m); % # of different signals (assumed to be in the same length(
M = zeros(W);

for k = 1:K
    
    for ii = 1:W+L-1
        
        if ii >= W-L && ii<=W
             M = M + m(k)*[zeros(W-ii,1); x(:,k); zeros(ii-L,1)]*[zeros(W-ii,1); x(:,k); zeros(ii-L,1)]';
        elseif ii < W-L
            M = M + m(k)*[zeros(W-ii,1); x(1:ii,k)]*[zeros(W-ii,1); x(1:ii,k)]';
        else
            M = M + m(k)*[x(ii-W+1:end,k); zeros(ii-L,1)]*[x(ii-W+1:end,k); zeros(ii-L,1)]';
        end
        
    end
end

M = M/N + sigma^2*eye(W);

end