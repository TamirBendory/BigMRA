function A = compute_A2(y_mat)

compute_ps = @(z) abs(fft(z)).^2;
W  = size(y_mat,1);
A = zeros(2*W-1,1);
N = size(y_mat,2);

for i = 1:N   
    
    a = ifft(compute_ps([y_mat(:,i) ; zeros(W-1,1)]));
    A =  A + a; %(1:W);
    
end
