function ps = psx(x,W,sigma,m,N)

% Computing the power spectrum
% x - signal
% W - window's length
% sigma - noise level
% m - # signals' repetitions
% N - measurement length

compute_ps = @(z) abs(fft(z)).^2;

L = length(x); % signal's length
K = length(m); % # of different signals (assumed to be in the same length(
%ps = zeros(W,1);
M = repmat(m',W,1);
ps = (W-L+1)*M.*compute_ps([x; zeros(W-L,K)]);

for ii = 0:L-2
    ps = ps + M.*compute_ps([x(1:ii+1,:) ;  zeros(W-ii-1,K)]);
end

for ii = 1:L-1
    ps = ps + M.*compute_ps([x(ii+1:end,:); zeros(W-L+ii,K)]);
end

ps = sum(ps,2)/N + sigma^2*W;

end



%
%
% for k = 1:K
%
%     ps = ps + m(k)*(W-L+1)*compute_ps([x(:,k); zeros(W-L,1)]);
%
%     for ii = 0:L-2
%         ps = ps + m(k)*compute_ps([x(1:ii+1,k) ;  zeros(W-ii-1,1)]);
%     end
%
%     for ii = 1:L-1
%         ps = ps + m(k)*compute_ps([x(ii+1:end,k); zeros(W-L+ii,1)]);
%     end
%
% end
