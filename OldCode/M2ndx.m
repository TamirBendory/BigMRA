function M = M2ndx(x,W,sigma,m,N,ac_flag)

% x - signal
% W - window's length
% sigma - noise level
% m - # signals' repetitions
% N - measurement length

L = length(x); % signal's length
K = length(m); % # of different signals (assumed to be of the same length)

if ac_flag
    M = zeros(W,1); M(1) = M(1) + sigma^2*W;
else
    M =  sigma^2*eye(W);
end

for k = 1:K
    
    if ac_flag
        
        for ii = W+1:W+L-1                        
                M = M + m(k)*x(ii-W+1,k)*[x(ii-W+1:end,k); zeros(ii-L,1)];       
        end
        
    else
        
        for ii = 1:W+L-1
        %for ii = W-L+1:W+L-1
            
            if ii >= W-L && ii<=W
                M = M + m(k)*[zeros(W-ii,1); x(:,k); zeros(ii-L,1)]*[zeros(W-ii,1); x(:,k); zeros(ii-L,1)]';
            elseif ii < W-L
                M = M + m(k)*[zeros(W-ii,1); x(1:ii,k)]*[zeros(W-ii,1); x(1:ii,k)]';
            else
                M = M + m(k)*[x(ii-W+1:end,k); zeros(ii-L,1)]*[x(ii-W+1:end,k); zeros(ii-L,1)]';
            end
            
        end
        
    end
end


M = M/N;

end