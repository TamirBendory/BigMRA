function [X_out, discrepancy_norm,err] = RRR_with_err(Y,L,th,X_init,X_true,hgram)

% The RRR algorithm estimates a WXW image X, supported on LXL pixels, from
% the measured Fourier magnitudes Y
% Note that Y = abs(fft2(X)) without the square
% X_true - the underlying signal

W = size(Y,1);
assert(W == size(Y,2),'The function supports only square images');
X_out = X_init(1:L,1:L);
% P1 - support projection
%Mask = ones(W);
%Mask(L+1:W,:) = zeros(W-L,W);
%Mask(1:L,L+1:W) = zeros(L,W-L);
load ('Mask'); % cheating by plugging exact mask
P1 = @(X) Mask.*X;

% P2 - Fourier magnitude projection
P2 = @(Z) Y.*sign(Z);

max_iter = 1e5; 
beta = 1;

% initial guess
if ~exist('X_init','var')
    X_init = zeros(W); X_init(1:L,1:L) = rand(L);
end

discrepancy_norm = zeros(max_iter,1);
err = discrepancy_norm;
X = X_init;

% iterations
for k = 1:max_iter
   
    X1 = P1(X); % support projection
    X1(X1<-0.0319) = -0.0319;
    X1(X1>1.0019) = 1.0019;
%    X1 = histeq(X1); %reshape(histeq(X1(:),hgram),size(X,1),size(X,2));
    err1 = norm(X1(1:L,1:L) - X_true,'fro')/norm(X_true(:));
    err2 = norm(rot90(X1(1:L,1:L),2) - X_true,'fro')/norm(X_true(:));
    err(k) = min(err1,err2);
    if and(err(k) < min(err(1:k-1)),k > 1) 
        X_out = X1(1:L,1:L); 
    end
    X2 = real(ifft2(P2(fft2(2*X1 - X)))); % Fourier magnitude projection
    discrepancy = X2 - X1;
    X = X + beta*discrepancy;% signal update
    %Xf = fft2(X);
    %Xf_phase = angle(Xf);
    %Xf_phase = histeq(angle(fft2(X_true_pad))); %reshape(histeq(Xf_phase,hgram),[size(X,1),size(X,2)]);
%    Xf = abs(Xf).*exp(1i*Xf_phase);
 %   X = real(ifft2(Xf));
    discrepancy_norm(k) = norm(discrepancy,'fro')/norm(X(:));    
    %X = reshape(histeq(X(:),hgram),size(X,1),size(X,2));
      
    if mod(k,1000) == 0
        fprintf('iter = %d, discrepancy = %.4g, err = %.4g\n',k,discrepancy_norm(k),err(k));
    end    
    if discrepancy_norm(k)<th
        fprintf('RRR last iteration = %d\n',k);
        break;
    end
end

X = P1(X);
X(X<0) = -0.0319;
X(X>1) = 1.0019;
%X = reshape(histeq(X(:),hgram),size(X,1),size(X,2));
err1 = norm(X(1:L,1:L) - X_true,'fro')/norm(X_true(:));
err2 = norm(rot90(X(1:L,1:L),2) - X_true,'fro')/norm(X_true(:));
err(k) = min(err1,err2);
if err(k) < min(err(1:k-1)) 
        X_out = X(1:L,1:L); 
end

if k == max_iter
    fprintf('RRR last iteration = %d\n',max_iter);
end

discrepancy_norm = discrepancy_norm(1:k);
err = err(1:k);

end

