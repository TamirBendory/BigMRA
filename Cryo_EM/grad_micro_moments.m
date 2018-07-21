function [g, store] = grad_micro_moments(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda, store)

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
    G_m2 = store.G_m2;
    G_m1 = store.G_m1;
    m_harms = store.m_harms;
end

if ~isfield(store, 'J')
    [G1, G2, G3] = bispectrum_grad_from_harmonics_par_noKloop...
        (a_lms, L, B, W, B_lists.L1, B_lists.L2, B_lists.L3, x_lists.L, x_lists.m, x_lists.s, B_lists.blk_id);
    
    sign_factor = (-1).^(x_lists.L(x_lists.n)+x_lists.m(x_lists.n));
    
    J = G1 + G2;
    J = bsxfun(@times, J(x_lists.n,:), sign_factor) + G3(x_lists.p, :);
    
    J = [lambda(1)*J, 2*lambda(2)*G_m2(x_lists.p,:)]; % add m2 derivative
    
    J(x_lists.m(x_lists.p) > 0,:) = 2*J(x_lists.m(x_lists.p) > 0, :); % treat m = 0 case
    
    J = [J, lambda(3)*G_m1(x_lists.p)]; % add m1 derivative
    J = gamma*[real(J); imag(J)].'; % separate real and imaginary parts
    
%     J = [J, m_harms]; % add gamma derivative
    store.J = J;
else
    J = store.J;
end

if ~isfield(store, 'g')
    g = 2*J.'*F;
    store.g = g;
else
   g = store.g;
end