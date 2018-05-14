function [X, discrepancy_norm] = RRR(Y,L,th,X_init,dc)

% The RRR algorithm estimates a WXW image X, supported on LXL pixels, from
% the measured Fourier magnitudes Y
% Note that Y = abs(fft2(X)) without the square

W = size(Y,1);
assert(W == size(Y,2),'The function supports only square images');

% P1 - support projection
%Mask = ones(W);
%Mask(L+1:W,:) = zeros(W-L,W);
%Mask(1:L,L+1:W) = zeros(L,W-L);
load ('Mask'); % cheating by plugging exact mask
P1 = @(X) Mask.*X;

% P2 - Fourier magnitude projection
P2 = @(Z) Y.*sign(Z);

max_iter = 1e8; %10^6;
beta = .5;

if exist('dc','var')
    dc_flag = 1;
else
    dc_flag = 0;
end

% initial guess
if ~exist('X_init','var')
X_init = zeros(W); X_init(1:L,1:L) = rand(L);
end
if dc_flag
X(1:L,1:L) = X(1:L,1:L)*dc/sum(X(:));    
end
discrepancy_norm = zeros(max_iter,1);
X = X_init;
% iterations
for k = 1:max_iter
    X1 = P1(X); % support projection
    %X1(X1<0) = 0; % for the case of image in [0,1]; should be removed otherwise
    X1(X1<-0.0319) = -0.0319; 
    %X1(X1>1) = 1;
    X1(X1>1.0019) = 1.0019; 
    X2 = real(ifft2(P2(fft2(2*X1 - X)))); % Fourier magnitude projection
    discrepancy = X2 - X1;
    X = X + beta*discrepancy;% signal update
    if dc_flag
    X(1:L,1:L) = X(1:L,1:L)*dc/sum(X(:));    
    end
    discrepancy_norm(k) = norm(discrepancy,'fro')/norm(X(:));
    
    if mod(k,1000) == 0
       fprintf('iter = %d, discrepancy = %.4g\n',k,discrepancy_norm(k));  
    end
    
    if discrepancy_norm(k)<th
        fprintf('RRR last iteration = %d\n',k);
        break;
    end
end

X = P1(X);
X(X<0) = 0;
X(X>1) = 1;

if dc_flag
X(1:L,1:L) = X(1:L,1:L)/sum(X(:))*dc;
end

if k == max_iter
    fprintf('RRR last iteration = %d\n',max_iter);
end

discrepancy_norm = discrepancy_norm(1:k);

end

