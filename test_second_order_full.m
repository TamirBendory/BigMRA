%% Forward model
L = 5;
W = 2*L;
x_true = randn(L, 1);
x_obs = [zeros(W-1, 1) ; x_true ; zeros(W-1, 1)];

M = zeros(W, W);
count = length(x_obs) - (W-1);
for k = 1 : count
    y = x_obs(k + (0:(W-1)));
    M = M + fft(y)*fft(y)';
end
M = M / count;

%% Inverse problem

params = struct();
params.L = L;
params.W = W;
params.M = M;

f = @test_second_order_full_cost;

manifold = euclideanfactory(L, 1);

% Needed to add this in reverse mode (it solved the problem I had then)
manifold.egrad2rgrad = @(x, g) real(g);

problem = manoptAD(manifold, f, params);

% Run it a few times from different random inits: you get a global opt
% pretty much every time, and at some point you get the right signal. Since
% minimizers appear to be isolated (see below), this suggests there are a
% finite number of global opts, and one of them is right. Just need a
% criterion to pick the right one.. (and a way to list them.) Count?
%
% The observation holds even with W = L (original test was W = 2L.)
x_est = trustregions(problem);

% hessianspectrum(problem, x_est) % -- suggests minimizers are isolated.

norm(x_true - x_est) / norm(x_true)
norm(x_true + x_est) / norm(x_true)

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