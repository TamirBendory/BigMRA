function [A, A_comp] = compute_A3(y_mat,L)

% I am taking L only for comparison

compute_bs = @(z) (fft(z)*(fft(z)')) .* circulant(fft(z));
W  = size(y_mat,1);
A = zeros(2*W);
N = size(y_mat,2);

for i = 1:N   
    
    a = real(ifft2(compute_bs([y_mat(:,i) ; zeros(W,1)])));
    %a = a(1:W,1:W);
    a = fftshift(a);
    A =  A + a; %(1:W,1:W);
   
  %  A =  A + a((W+1)/2:(3*W+1)/2-1,(W+1)/2:(3*W+1)/2-1);
    
end

%A = A + flipud(A);
%A = A + fliplr(A);
A_comp = A;
%Q = (linspace(ceil(L/2),L,L))';
%Q = [flipud(Q) ; ones(W-L,1)];
%A_comp = A/N./Q/2;