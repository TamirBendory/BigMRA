%% Forward model
L = 5;
W = 10;
x_true = randn(L, 1);
x_obs = [zeros(W-1, 1) ; x_true ; zeros(W-1, 1)];

M = zeros(W, 1);
count = length(x_obs) - (W-1);
for k = 1 : count
    y = x_obs(k + (0:(W-1)));
    M = M + abs(fft(y)).^2;
end
M = M / count;

%% Inverse problem

params = struct();
params.L = L;
params.W = W;
params.M = M;

f = @test_second_order_cost;

manifold = euclideanfactory(L, 1);

problem = manoptAD(manifold, f, params);

x_est = trustregions(problem);

norm(x_true - x_est) / norm(x_true)
norm(x_true + x_est) / norm(x_true)

%%
return;

%% Full second order
L = 5;
W = L;
x_true = randn(L, 1);
x_obs = [zeros(W-1, 1) ; x_true ; zeros(W-1, 1)];

M = zeros(W, W);
count = length(x_obs) - (W-1);
for k = 1 : count
    y = x_obs(k + (0:(W-1)));
    M = M + y*y';
end
M = M / count;
imagesc(M); axis equal; colorbar;