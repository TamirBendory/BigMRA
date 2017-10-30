function x = my_pinv(hy,hs,eps)

x = hy./hs;
ind = find (abs(hs)<eps);
x(ind) = 0;