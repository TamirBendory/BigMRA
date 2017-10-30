clear; close all; clc;

L = 10; % length of signal
x = randn(L,1); 

N = 200000; 
k = 20000;
sigma = 0;
window_size = L; 

[y ind] = gen_data(x,N,k,L,sigma); 
s = zeros(N,1); s(ind) = 1;
%% rearraging the matrix data

M = N/window_size;
y_matrix = zeros(M,window_size);
y_stretch = [y ; y(1:window_size-1)];
for i = 1:M
    y_matrix(i,:) = y_stretch(i:i+window_size-1);
end
clear y_stretch;

%% invariants

[mean_est, P_est, B_est] = invariants_from_data(y_matrix', sigma);
[z, problem] = phases_from_bispectrum_real(B_est, sign(mean_est), randn(window_size,1));

x_est = real(ifft(sqrt(P_est).*z));
xref = zeros(window_size,1); xref(1:L) = x;

[x_aligned,~, ~] = align_to_reference(x_est, xref);

% thresholding

err = norm(x_aligned(1:L) - x)/norm(x);

figure; hold on;  
stem(xref); stem(x_aligned,'*r'); 
legend('signal','recovered');