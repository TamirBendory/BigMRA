function [f, store] = cost_micro_moments(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda, store)

% gamma = x(end);
gamma = 1;
% a_lms = x(1:(end-1)/2) + 1i*x((end-1)/2+1:end-1);
a_lms = x(1:end/2) + 1i*x(end/2+1:end);
a_lms = vec_to_cell_vol_coeffs(a_lms, x_lists.L(x_lists.p));
a_vec = cellfun(@(x)x(:), a_lms, 'UniformOutput', 0);
a_vec = vertcat(a_vec{:});

if ~isfield(store, 'F')
    m3_harms = bispectrum_from_harmonics_par_noKloop(a_lms, L, B, W, B_lists.L1, B_lists.L2, B_lists.L3, x_lists.s, B_lists.blk_id);
    dm3 = gamma*m3_harms - m3_micro;
    
    [m2_harms, G_m2] = power_spectrum_from_harmonics(a_vec, M2_quants.j_l, M2_quants.R_0n, M2_quants.w, L);
    dm2 = gamma*m2_harms - m2_micro;
    
    [m1_harms, G_m1] = mean_from_harmonics(a_vec, M2_quants.j_0, x_lists.s(1), L);
    dm1 = gamma*m1_harms - m1_micro;
    
    F = real([lambda(1)*dm3; lambda(2)*dm2; lambda(3)*dm1]);
    store.m_harms = real([lambda(1)*m3_harms; lambda(2)*m2_harms; lambda(3)*m1_harms]);
    store.F = F;
    store.G_m2 = G_m2;
    store.G_m1 = G_m1;
else
    F = store.F;
end

if ~isfield(store, 'f')
    f = norm(F)^2;
    store.f = f;
else
    f = store.f;
end