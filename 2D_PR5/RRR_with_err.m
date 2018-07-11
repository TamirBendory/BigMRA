function [X_out, discrepancy_norm,err,err1,err2] = RRR_with_err(Y,L,th,X_init,X_true)

% The RRR algorithm estimates a WXW image X, supported on LXL pixels, from
% the measured Fourier magnitudes Y
% Note that Y = abs(fft2(X)) without the square
% X_true - the underlying signal
% th - threshold to halt the RRR algorithm

W = size(Y,1);
assert(W == size(Y,2),'The function supports only square images');
X_out = X_init(1:L,1:L);
% P1 - support projection
Mask = ones(W);
Mask(L+1:W,:) = zeros(W-L,W);
Mask(1:L,L+1:W) = zeros(L,W-L);
P1 = @(X) Mask.*X;

% P2 - Fourier magnitude projection
P2 = @(Z) Y.*sign(Z);

max_iter = 1000; 
beta = 1;

% initial guess
if ~exist('X_init','var')
    X_init = zeros(W); X_init(1:L,1:L) = rand(L);
end

discrepancy_norm = zeros(max_iter,1);
err = discrepancy_norm;
err1 = err;
err2 = err;
X = X_init;

% iterations
for k = 1:max_iter
    X1 = P1(X); % support projection
    X1(X1<-1) = -1;
    X1(X1>1) = 1;
    err1(k) = norm(X1(1:L,1:L) - X_true,'fro')/norm(X_true(:));
    err2(k) = norm(rot90(X1(1:L,1:L),2) - X_true,'fro')/norm(X_true(:));
    err(k) = min(err1(k),err2(k));
    if and(err(k) < min(err(1:k-1)),k > 1) 
        X_out = X1(1:L,1:L); 
    end
    X2 = real(ifft2(P2(fft2(2*X1 - X)))); % Fourier magnitude projection
    discrepancy = X2 - X1;
    X = X + beta*discrepancy;% signal update
    discrepancy_norm(k) = norm(discrepancy,'fro')/norm(X(:));          
%    if mod(k,1000) == 0
%        fprintf('iter = %d, discrepancy = %.4g, err = %.4g\n',k,discrepancy_norm(k),err(k));
%    end    
%    if discrepancy_norm(k)<th
%        fprintf('RRR last iteration = %d\n',k);
%        break;
%    end
end

X = P1(X);
X(X<-1) = 1;
X(X>1) = 1;
err1(k) = norm(X(1:L,1:L) - X_true,'fro')/norm(X_true(:));
err2(k) = norm(rot90(X(1:L,1:L),2) - X_true,'fro')/norm(X_true(:));
err(k) = min(err1(k),err2(k));
if err(k) < min(err(1:k-1)) 
       X_out = X(1:L,1:L); 
end

%if k == max_iter
%    fprintf('RRR last iteration = %d\n',max_iter);
%end

discrepancy_norm = discrepancy_norm(1:k);
err = err(1:k);

end

