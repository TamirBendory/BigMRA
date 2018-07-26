% Test matrix-matrix multiplication vs. matrix-column mult.

n = 1e3; m = 1e5;
M = randn(n); % matrix to be applied
rng(1)
tic,
A = zeros(n, m, 'double');
for ii = 1:m
    v = randn(n, 1);
    A(:, ii) = M*v;
end
time_MatCol = toc

rng(1)
tic,
B = zeros(n, m, 'double');
for ii = 1:m
    v = randn(n, 1);
    B(:, ii) = v;
end 
A2 = M*B;
time_MatMat = toc

norm(A(:) - A2(:))/norm(A(:))