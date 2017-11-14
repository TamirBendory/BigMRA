function x_aligned = auto_alignment(x_est,L,flag_plot,x)

% input: x_est - the estimated signal, which is rotated and padded by W-L zeros
%        L - the length of the output
%        flag_plot - plots the alignment process if 1 
%        x - the true signal - only for comparison

% output: the cropped aligned signal

% computing the energy at each window
W = length(x_est);
eng = zeros(W,1);
for i = 0:W-1
    eng(i+1) = norm(x_est(i+1:mod(i+L-1,W)+1));
end

% alignment
[~,est_shift] = max(eng);
x_aligned = circshift(x_est,-est_shift+1);


if flag_plot
figure; 
subplot(211);hold on;
plot(x_aligned); plot(x_est); plot([x ; zeros(W-L,1)]);
legend('aligned estimation','estimated signal','true signal');
subplot(212); plot(eng); title('the alignment map');
end

% cropping
x_aligned = x_aligned(1:L);
