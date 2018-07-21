function [f, store] = cost_pwr_spectrum(a_lms, x_lists, L, M2_quants, m2_micro, gamma, store)

a_lms = a_lms(1:length(a_lms)/2) + 1i*a_lms(length(a_lms)/2+1:end);
a_lms = vec_to_cell_vol_coeffs(a_lms, x_lists.L(x_lists.p));
a_vec = cellfun(@(x)x(:), a_lms, 'UniformOutput', 0);
a_vec = vertcat(a_vec{:});

if ~isfield(store, 'dm2')
    [m2_harms, G_m2, H_m2] = power_spectrum_from_harmonics(a_vec, M2_quants.j_l, M2_quants.R_0n, M2_quants.w, L);
    dm2 = m2_harms - m2_micro./gamma;
    store.dm2 = dm2;
    store.G_m2 = G_m2;
    store.H_m2 = H_m2;
else
    dm2 = store.dm2;
end

if ~isfield(store, 'f')
    f = norm(dm2)^2;
    store.f = f;
else
    f = store.f;
end