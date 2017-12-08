% launch test_one_signal_in_micrograph
x_zp = [x_true ; zeros(W-L, 1)];
% [ifft(abs(fft(x_zp)).^2) , [M2(1:(W+1)/2) ; flipud(M2(2:(W+1)/2))]]

% This is the power spectrum of x zero padded to length W
[abs(fft(x_zp)).^2, fft([M2(1:(W+1)/2) ; flipud(M2(2:(W+1)/2))])]

% This works if M3 computed from -(W-1)/2 : (W-1)/2 with W odd
B_x_zp = (fft(x_zp)*fft(x_zp)') .* circulant_AD(fft(x_zp));
subplot(1,2,1);
imagesc(ifft2(B_x_zp));
subplot(1,2,2);
imagesc(rot90(rot90(fftshift(M3))));

% B_x_zp = (fft(x_zp)*fft(x_zp)') .* circulant_AD(fft(x_zp));
% subplot(1,2,1);
% imagesc(imag(B_x_zp));
% subplot(1,2,2);
% imagesc(imag(fft2(rot90(rot90(fftshift(M3))))));


% subplot(2,2,3);
% imagesc((real(B_x_zp)));
% subplot(2,2,4);
% imagesc(real(fft2(M3')));
