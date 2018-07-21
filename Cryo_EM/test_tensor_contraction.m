function [t1, t2] = test_tensor_contraction(L)
% compute sum_{i,j,k} A_{i,j,k} * x_i * x_j * x_k

A = randn(L, L, L);
x = randn(L,1);

% Method 1:
tic
tmp = reshape(A, L^2, L)*x;
tmp = reshape(tmp, L, L)*x;
tmp = tmp.'*x;
t1 = toc;

% Method 2:
tic
x = reshape(x*x.', L^2, 1)*x.';
% x = kron(kron(x, x), x);
tmp2 = A(:).'*x(:);
t2 = toc;

assert(abs(tmp - tmp2)<1e-10)