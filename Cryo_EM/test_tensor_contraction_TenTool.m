% Test tensor contraction using Sandia Lab's Tensor Toolbox
addpath('../tensor_toolbox')
addpath('../tensor_lab')
addpath('../tprod')

X_mat = randn(50,30,40,60,70);

A = randn(70,50);
B = randn(60,40);
C = randn(40,40);

% Naive contraction:
tic
Res_naive = reshape(X_mat, [], size(A,1))*A; % 5*3*4*6 x 5
Res_naive = reshape(Res_naive.', [], size(B,1))*B;
Res_naive = reshape(Res_naive.', [], size(C,1))*C;
Res_naive = reshape(Res_naive.', size(C,2), size(B,2), size(A,2), size(X_mat,1), size(X_mat,2));
Res_naive = permute(Res_naive, [4,5,1,2,3]);
time_naive = toc;

% Tensor Toolbox:
X = tensor(X_mat); % 5 x 3 x 4 x 2

tic
Res_TenTool = ttm(X, {C.',B.',A.'}, [3,4,5]); 
time_TenTool = toc;

% Tensor Lab:
tic
Res_TenLab = tmprod(X_mat, {C.', B.', A.'}, 3:5);
time_TenLab = toc;

% Tprod:
tic
Res_tprod = tprod(X_mat, [1,2,3,4,-1], A, [-1, 5]);
Res_tprod = tprod(Res_tprod, [1,2,3,-1,5], B, [-1, 4]);
Res_tprod = tprod(Res_tprod, [1,2,-1,4,5], C, [-1, 3]);
time_tprod = toc;

assert(norm(Res_naive(:)-Res_TenTool.data(:)) < 1e-10)
assert(norm(Res_naive(:)-Res_TenLab(:)) < 1e-10)