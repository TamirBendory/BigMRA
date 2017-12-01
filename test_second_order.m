%% Forward model
L = 5;
W = 2*L-1;
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
% Anytime complex numbers might appear you need this for AD ..
manifold.egrad2rgrad = @(x, g) real(g);

problem = manoptAD(manifold, f, params);

x_est = trustregions(problem);

norm(x_true - x_est) / norm(x_true)
norm(x_true + x_est) / norm(x_true)

% Same observation as for the full second order test: the true signal is an
% isolated global optimizer, at least for W = 2L, and we can get it (by
% chance.)

%%
return;

%% Full second order
L = 10;
W = 2*L;
x_true = randn(L, 1);
x_obs = [zeros(W-1, 1) ; x_true ; zeros(W-1, 1)];

order = 2;
M = zeros(W^order, 1);
count = length(x_obs) - (W-1);
for k = 1 : count
    y = fft(x_obs(k + (0:(W-1))));
    m = y;
    for q = 2 : order
        m = kron(m, conj(y)); % though we should also look at conjugates...
    end
    M = M + m;
end
M = M / count;
numel(uniquetol(abs([real(M) ; imag(M)]), 1e-10))


subplot(1,2,1);
imagesc(abs(real(reshape(M, [W, W]))));
axis equal; axis tight; colorbar;
subplot(1,2,2);
imagesc(abs(imag(reshape(M, [W, W]))));
axis equal; axis tight; colorbar;