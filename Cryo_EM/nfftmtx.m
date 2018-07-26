function F = nfftmtx(N, x)
% NFFT matrix for an N x N image with nodes x. 

if mod(N,2) == 0
    k = -N/2:N/2-1;
else
    k = -(N-1)/2:(N-1)/2;
end
[K1,K2]=meshgrid(k,k);
K1=K1(:);
K2=K2(:);

F = exp(-1i*2*pi*(x(:,1)*K1' + x(:, 2)*K2'));