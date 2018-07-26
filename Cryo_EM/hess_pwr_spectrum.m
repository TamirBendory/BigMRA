function [Hu, store] = hess_pwr_spectrum(a_lms, x_lists, L, M2_quants, m2_micro, gamma, u, store)

a_lms = a_lms(1:length(a_lms)/2) + 1i*a_lms(length(a_lms)/2+1:end);
a_lms = vec_to_cell_vol_coeffs(a_lms, x_lists.L(x_lists.p));
a_vec = cellfun(@(x)x(:), a_lms, 'UniformOutput', 0);
a_vec = vertcat(a_vec{:});

if ~isfield(store, 'dm2')
    [m2_harms, G_m2, H_m2] = power_spectrum_from_harmonics(a_vec, M2_quants.j_l, M2_quants.R_0n, M2_quants.w, L);
    dm2 = real(m2_harms - m2_micro./gamma);
    store.dm2 = dm2;
    store.G_m2 = G_m2;
    store.H_m2 = H_m2;
else
    dm2 = store.dm2;
    G_m2 = store.H_m2;
    H_m2 = store.H_m2;
end

if ~isfield(store, 'J')
    J = 2*G_m2(x_lists.p,:);
    J(x_lists.m(x_lists.p) > 0,:) = 2*J(x_lists.m(x_lists.p) > 0, :);
    
    J = [real(J); imag(J)].';
    store.J = J;
else
    J = store.J;
end

if ~isfield(store, 'H')
    H_m2 = cellfun(@(x) x(x_lists.p, x_lists.p), H_m2, 'UniformOutput', 0);
    for ii = 1:length(H_m2)
        H_m2{ii}(:, x_lists.m(x_lists.p) > 0) = 2*H_m2{ii}(:, x_lists.m(x_lists.p) > 0);
    end
    m_zero_inds = x_lists.m(x_lists.p)==0;
    sign_factor = (-1).^(x_lists.L(x_lists.p) + x_lists.m(x_lists.p));
    H_m2_real = H_m2; H_m2_imag = H_m2;
    for ii = 1:length(H_m2)
        H_m2_real{ii}(m_zero_inds, m_zero_inds) = (H_m2_real{ii}(m_zero_inds, m_zero_inds) + ...
            bsxfun(@times, H_m2_real{ii}(m_zero_inds, m_zero_inds), sign_factor(m_zero_inds)))/2;
        H_m2_imag{ii}(m_zero_inds, m_zero_inds) = (H_m2_imag{ii}(m_zero_inds, m_zero_inds) - ...
            bsxfun(@times, H_m2_imag{ii}(m_zero_inds, m_zero_inds), sign_factor(m_zero_inds)))/2;
    end
    H = cellfun(@(x,y) 2*[real(x), zeros(size(x)); zeros(size(x)), real(y)], H_m2_real, H_m2_imag, 'UniformOutput', 0);
else
    H = store.H;
end

h = 0;
for ii = 1:length(H)
    h = h + dm2(ii)*H{ii};
end
h = h + J.'*J;
h = 2*h;
store.h = h;

% u = cell_to_vec_vol_coeffs(vec_to_cell_vol_coeffs(u(1:end/2)+1i*u(end/2+1:end), x_lists.L(x_lists.p)));
% u = [real(u);imag(u)];
Hu = h*u;