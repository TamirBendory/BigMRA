% testing the RRR

W = 13;
L = 5;
X = zeros(W);
X(1:L,1:L) = rand(L);
Y = abs(fft2(X));

Xout = RRR(Y,L);
Xc = X(1:L,1:L); Xout = Xout(1:L,1:L);
figure; imagesc([Xc,Xout]);
err = norm(Xc(:) - Xout(:))/norm(Xc(:))