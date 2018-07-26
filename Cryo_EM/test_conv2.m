% Test FFT-based 2D convolution against brute-force sliding window

L = 31
M = 8e3
I = randn(M, M); % "micrograph"
W = randn(2*L-1, 2*L-1); % sliding window

% Fast implementation:
addpath('../CONVNFFT_Folder')
% C = convnfft(I, W, 'same');
tic
C = xcorr2_fft(I, W);
C = C(L:end-L+1, L:end-L+1);
time_fast = toc

% MEX implementation:
opts.GPU = false;
opts.Power2Flag = true;
tic
C3 = convnfft(I, W(end:-1:1, end:-1:1), 'same', 1:2, opts);
time_convnfft = toc

% Brute force:
C2 = zeros(size(I));
I = padarray(I, [L-1, L-1]);
tic
for col = L:M+L-1
    for row = L:M+L-1
        curr_patch = I(row-L+1:row+L-1, col-L+1:col+L-1);
        C2(row-L+1, col-L+1) = dot(curr_patch(:), W(:));
    end
end
time_brute = toc

norm(C(:)-C2(:))/norm(C2(:))
