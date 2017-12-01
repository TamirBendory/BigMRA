function A = A2x(x,W,sigma,m,N)

% Computing the power spectrum
% x - signal
% W - window's length
% sigma - noise level
% m - # signals' repetitions
% N - measurement length

compute_ps = @(z) abs(fft(z)).^2;

L = length(x); % signal's length
K = length(m); % # of different signals (assumed to be in the same length(
M = repmat(m',W,1);
s = ifft(compute_ps([x; zeros(W-L,K) ; zeros(W,K)]));
A = (W-L+1)*M.*s(1:W,:);

for ii = 0:L-2
    
    s = ifft(compute_ps([x(1:ii+1,:) ;  zeros(W-ii-1,K) ; zeros(W,K)]));
    A = A + M.*s(1:W,:);
    
end

for ii = 1:L-1
    
    s = ifft(compute_ps([x(ii+1:end,:); zeros(W-L+ii,K) ; zeros(W,K)]));
    A = A + M.*s(1:W,:);
    
end
% 
% Q = (linspace(ceil(L/2),L,L))';
% Q = [flipud(Q) ; ones(W-L,1)];
% A = sum(A,2)/N./Q/2;
A(1) = A(1) + sigma^2;

end


