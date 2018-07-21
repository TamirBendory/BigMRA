function [a_lms, gamma] = recover_vol_coeffs_from_moments(m3_micro, m2_micro, m1_micro, maxL, L, B, W, B_lists, r_cut, a_init)

% Get vectorization indices for the coefficients:
s_lens = gen_s_list(maxL, r_cut, 1, floor(L/2));
[L_list, m_list, ~, p_inds, n_inds] = gen_vec_coeff_lists(maxL, s_lens);
x_lists.L = L_list;
x_lists.m = m_list;
x_lists.p = p_inds;
x_lists.n = n_inds;
x_lists.s = s_lens;

% Get bases and quadrature for 1st and 2nd moment computations:
q0 = sqrt(sum(B_lists.blk_id == 1));
[r,w] = lgwt(20*q0, 0, 1);
w = w.*r;
j_l = generate_spherical_bessel_basis(maxL, s_lens, 1/2, (1/2)*r);
j_0 = cell2mat(generate_spherical_bessel_basis(0, s_lens, 1/2, 0));
[R_0n, alpha_Nn_2D] = PSWF_radial_2D(0, q0-1, pi*(L-1), r); % generate 2D radial prolates
R_0n = bsxfun(@times, R_0n, 2./alpha_Nn_2D(:).');

M2_quants.w = w;
M2_quants.j_l = j_l;
M2_quants.j_0 = j_0;
M2_quants.R_0n = R_0n;

% Initialization:
if ~exist('a_init','var') || isempty(a_init)
    a_init = randn(2*length(p_inds), 1);
    a_init(end+1) = 0.5;
end

% Bounds for occupancy factor:
lb = zeros(2*length(p_inds)+1, 1, 'double');
lb(1:end-1) = -inf;
lb(end) = 0;

ub = zeros(2*length(p_inds)+1, 1, 'double');
ub(1:end-1) = inf;
ub(end) = 1;

% Weights for least-squares (lambda(1) -> 3rd moment, lambda(2) -> 2nd
% moment, lambda(3) -> 1st moment)
% lambda(1) = 1/(norm(m3_micro)); 
% lambda(2) = 1/(norm(m2_micro));
% lambda(3) = 1/abs(m1_micro);

lambda(1) = 1/(norm(m3_micro))/numel(m3_micro);
lambda(2) = 1/(norm(m2_micro))/numel(m2_micro);
lambda(3) = 1/(norm(m1_micro))/numel(m1_micro);

costgrad = @(x) costgrad_lsqnonlin(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda);
options = optimoptions('lsqnonlin','Jacobian','on','DerivativeCheck','off','Display','iter','TolX',1e-8,'TolFun',1e-8,'MaxIter',1e4);
tic, x = lsqnonlin(costgrad, a_init, lb, ub, options); time_TR = toc
% tic, x = lsqnonlin(costgrad, a_init, [], [], options); time_TR = toc

% problem.M = euclideanfactory(2*length(p_inds)+1, 1);
% 
% problem.costgrad = @(x,store) costgrad_micro_moments(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda, store);
% problem.hess = @(x, u, store) hess_micro_moments(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda, u, store);
% % checkgradient(problem)
% % checkhessian(problem)
% 
% % options.maxiter = 1e4;
% opts_manopt.tolgradnorm = 1e-5;
% % options.tolcost = 1e-12;
% % % options.statsfun = @statsfun_bispect;
% tic, x = trustregions(problem, a_init, opts_manopt); time_manopt = toc

gamma = x(end);
a_lms = x(1:(end-1)/2) + 1i*x((end-1)/2+1:end-1);
a_lms = vec_to_cell_vol_coeffs(a_lms, L_list(m_list>=0));
