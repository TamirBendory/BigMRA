function [Hu, store] = hess_micro_moments_withHess(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda, u, store)

if ~isfield(store, 'F')
    [~, store] = cost_micro_moments_withHess(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda, store);
end
H1 = store.H1;
H2 = store.H2;
H3 = store.H3;
H_m2 = store.H_m2;
F = store.F;

if ~isfield(store, 'J')
    [~, store] = grad_micro_moments_withHess(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda, store);
end
J = store.J;
sign_factor = store.signs;

if ~isfield(store, 'h')
    Hc = cellfun(@(x,y,z) (sign_factor*sign_factor.').*(z(x_lists.n,x_lists.n) + z(x_lists.n,x_lists.n).') ...
        + bsxfun(@times, x(x_lists.p, x_lists.n), sign_factor.') + bsxfun(@times, x(x_lists.p, x_lists.n).', sign_factor)...
        + bsxfun(@times, y(x_lists.p, x_lists.n), sign_factor.') + bsxfun(@times, y(x_lists.p, x_lists.n).', sign_factor)...
        , H1, H2, H3, 'UniformOutput', 0);
    
    Hr = cellfun(@(x,y,z) bsxfun(@times, z(x_lists.p, x_lists.n), sign_factor.') + bsxfun(@times, z(x_lists.n, x_lists.p), sign_factor).' ...
        + (sign_factor*sign_factor.').*x(x_lists.n,x_lists.n) + x(x_lists.p, x_lists.p).'...
        + (sign_factor*sign_factor.').*y(x_lists.n,x_lists.n) + y(x_lists.p, x_lists.p).'...
        , H1, H2, H3, 'UniformOutput', 0);
    
    H_plus = cellfun(@plus, Hc, Hr, 'UniformOutput', 0);
    H_minus = cellfun(@minus, Hc, Hr, 'UniformOutput', 0);
    for ii = 1:length(H_plus)
        H_plus{ii}(:, x_lists.m(x_lists.p) > 0) = 2*H_plus{ii}(:, x_lists.m(x_lists.p) > 0);
        H_plus{ii}(x_lists.m(x_lists.p) == 0, :) = (1/2)*H_plus{ii}(x_lists.m(x_lists.p) == 0, :);
        
        H_minus{ii}(:, x_lists.m(x_lists.p) > 0) = 2*H_minus{ii}(:, x_lists.m(x_lists.p) > 0);
        H_minus{ii}(x_lists.m(x_lists.p) == 0, :) = (1/2)*H_minus{ii}(x_lists.m(x_lists.p) == 0, :);
    end
    H = cellfun(@(x,y) [real(x), imag(x); imag(y), -real(y)], H_plus, H_minus, 'UniformOutput', 0);
    
    H_m2 = cellfun(@(x) x(x_lists.p, x_lists.p), H_m2, 'UniformOutput', 0);
    for ii = 1:length(H_m2)
        H_m2{ii}(:, x_lists.m(x_lists.p) > 0) = 2*H_m2{ii}(:, x_lists.m(x_lists.p) > 0);
    end
    H_m2 = cellfun(@(x) 2*lambda*[real(x), zeros(length(x_lists.p),'double');...
        zeros(length(x_lists.p), 'double'), real(x)], H_m2, 'UniformOutput', 0);
    H = [H(:);H_m2(:)];
    
    h = 0;
    for ii = 1:length(H)
        h = h + F(ii)*H{ii};
    end
    h = h + J.'*J;
    h = 2*h;
    store.h = h;
else
    h = store.h;
end
u = cell_to_vec_vol_coeffs(vec_to_cell_vol_coeffs(u(1:end/2)+1i*u(end/2+1:end), x_lists.L(x_lists.p)));
u = [real(u);imag(u)];
Hu = h*u;
