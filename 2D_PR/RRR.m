function X = RRR(Y,L)

% The RRR algorithm estimates a WXW image X, supported on LXL pixels, from
% the measured Fourier magnitudes Y 
% Note that Y = abs(fft2(X)) without the square

W = size(Y,1);
assert(W == size(Y,2),'The function supports only square images');

% The algorithm builds upon two building blocks
% P1 - support projection
Mask = ones(W); 
Mask(L+1:W,:) = zeros(W-L,W);
Mask(1:L,L+1:W) = zeros(L,W-L);    
P1 = @(X) Mask.*X; 
% P2 - Fourier magnitude projection
P2 = @(Z) Y.*sign(Z);

max_iter = 10^4;
beta = 0.5;
% stopping criterion
th = 10^-8;
% initial guess
X = zeros(W); X(1:L,1:L) = rand(L); 

% iterations
for k = 1:max_iter
   X1 = P1(X); % support projection
   X2 = real(ifft2(P2(fft2(2*X1 - X)))); % Fourier magnitude projection
   discrepancy = X2 - X1; 
   X= X + beta*discrepancy;% signal update
   if norm(discrepancy,'fro')<th
       fprintf('RRR last iteration = %d\n',k);
       break;
   end
end
X = P1(X);
if k == max_iter
       fprintf('RRR last iteration = %d\n',max_iter);    
end
X = X(1:L,1:L);
end

