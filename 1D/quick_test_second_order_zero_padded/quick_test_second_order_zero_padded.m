clear;
% See notes Jan 31, 2018
L = 3;
zp = @(x) [x ; zeros(L-1, 1)];
zp_adj = @(y) y(1:L);
x_true = randn(L, 1);
q = abs(fft(zp(x_true))).^2;
problem.M = euclideanfactory(L, 1); % but it's not quite what I wanted to try: I wanted to have W = 2L-1 here instead of L; though it shouldn't change anything because q indeed contains L distinct real numbers, and the optimization manages to match them perfectly, so it's not an optimization issue.
R = @(x) abs(fft(zp(x))).^2 - q;
problem.cost = @(x) .25/(2*L-1)*norm(R(x), 2)^2;
problem.grad = @(x) real(zp_adj(ifft(R(x).*fft(zp(x)))));
% checkgradient(problem);
x = trustregions(problem);
norm(x-x_true)/norm(x_true) % x is far from x_true with high probability
[ifft(abs(fft(zp(x_true))).^2), ifft(abs(fft(zp(x))).^2)] % but their moments match: second order is not enough
% this is a classical result in 1D phase retrieval
% I think there are supposed to be only a finite number of solutions, but
% it's an exponentially large number.