function bs = bsx(x,W,sigma,m,N,M1)

% x - signal
% W - window's length
% sigma - noise level
% m - # signals' repetitions
% N - measurement length
% M1 - the first moment (I am considering the real one; not the estimated)


%fftz = @(x) fft(z);
compute_bs = @(z) (fft(z)*(fft(z)')) .* circulant(fft(z));
L = length(x); % signal's length
K = length(m); % # of different signals (assumed to be in the same length(
%ps = zeros(W,1);
M = repmat(m',W,1);
bs = zeros(W);
A = eye(W); A(1,:) = A(1,:)+1; A(:,1) = A(:,1)+1;

for k = 1:K
    
    bs = bs + m(k)*(W-L+1)*compute_bs([x(:,k); zeros(W-L,1)]);
    
    for ii = 0:L-2
        bs = bs + m(k)*compute_bs([x(1:ii+1,k) ;  zeros(W-ii-1,1)]);
    end
    
    for ii = 1:L-1
        bs = bs + m(k)*compute_bs([x(ii+1:end,k); zeros(W-L+ii,1)]);
    end
    
end

bs = bs/N + sigma^2*W^2*A*M1;

end