function [f, store] = cost_micro_moments_withHess(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda, store)

a_lms = vec_to_cell_vol_coeffs(a_lms(1:length(a_lms)/2) + 1i*a_lms(length(a_lms)/2+1:end), x_lists.L(x_lists.p));
a_lms_vec = cellfun(@(x)x(:), a_lms, 'UniformOutput', 0);
a_lms_vec = vertcat(a_lms_vec{:});

if ~isfield(store, 'F')
    [H1, H2, H3] = bispectrum_hess_from_harmonics(a_lms, L, B, W, B_lists.L1, B_lists.L2, B_lists.L3, x_lists.L, x_lists.m, x_lists.s, B_lists.blk_id);

    [m2_harms, G_m2, H_m2] = power_spectrum_from_harmonics(a_lms_vec, M2_quants.j_l, M2_quants.R_0n, M2_quants.w, L);
    
    G1 = cellfun(@(x) x.'*a_lms_vec, H3, 'UniformOutput', 0);
    G1 = cat(2, G1{:});
    
    G2 = cellfun(@(x) x*a_lms_vec, H3, 'UniformOutput', 0);
    G2 = cat(2, G2{:});
    
    G3 = cellfun(@(x) x*a_lms_vec, H1, 'UniformOutput', 0);
    G3 = cat(2, G3{:});
    
    m3_harms = G1.'*a_lms_vec;

    dm3 = m3_harms - (1/gamma)*m3_micro;
    dm2 = m2_harms - (1/gamma)*m2_micro;
    F = [real(dm3); lambda*real(dm2)];
    
    f = norm(F)^2;
    
    store.H1 = H1;
    store.H2 = H2;
    store.H3 = H3;
    store.H_m2 = H_m2;
    store.G1 = G1;
    store.G2 = G2;
    store.G3 = G3;
    store.G_m2 = G_m2;
    store.F = F;
else
    f = norm(store.F)^2;
end


