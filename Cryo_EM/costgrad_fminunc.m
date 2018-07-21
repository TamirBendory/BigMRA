function [f, g, h] = costgrad_fminunc(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda)
% Calculate objective f
a_lms = a_lms(1:length(a_lms)/2) + 1i*a_lms(length(a_lms)/2+1:end);
a_lms = vec_to_cell_vol_coeffs(a_lms, x_lists.L(x_lists.p));

a_lms_vec = cellfun(@(x)x(:), a_lms, 'UniformOutput', 0);
a_lms_vec = vertcat(a_lms_vec{:});
[m2_harms, G_m2, H_m2] = power_spectrum_from_harmonics(a_lms_vec, M2_quants.j_l, M2_quants.R_0n, M2_quants.w, L);
dm2 = m2_harms - m2_micro./gamma;

switch nargout
    case 1 % only cost
        m3_harms = bispectrum_from_harmonics_par_noKloop(a_lms, L, B, W, B_lists.L1, B_lists.L2, B_lists.L3, x_lists.s, B_lists.blk_id);
        dm3 = m3_harms - (1/gamma)*m3_micro;
        
        F = [real(dm3); lambda*real(dm2)];
        f = norm(F)^2;
    case 2 % cost & grad
        [G1, G2, G3] = bispectrum_grad_from_harmonics_par_noKloop...
            (a_lms, L, B, W, B_lists.L1, B_lists.L2, B_lists.L3, x_lists.L, x_lists.m, x_lists.s, B_lists.blk_id);
        m3_harms = G1.'*a_lms_vec;
        
        dm3 = m3_harms - (1/gamma)*m3_micro;
        F = [real(dm3); lambda*real(dm2)];
        f = norm(F)^2;
        
        sign_factor = (-1).^(x_lists.L(x_lists.n)+x_lists.m(x_lists.n));
        
        J = G1 + G2;
        J = bsxfun(@times, J(x_lists.n,:), sign_factor) + G3(x_lists.p, :);
        J = [J, 2*lambda*G_m2(x_lists.p,:)];
        J(x_lists.m(x_lists.p) > 0,:) = 2*J(x_lists.m(x_lists.p) > 0, :);
        
        J = [real(J); imag(J)].';
        g = 2*J.'*F;
        
    case 3 % cost, grad, hess
        [H1, H2, H3] = bispectrum_hess_from_harmonics(a_lms, L, B, W, B_lists.L1, B_lists.L2, B_lists.L3, x_lists.L, x_lists.m, x_lists.s, B_lists.blk_id);
        
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
        
        sign_factor = (-1).^(x_lists.L(x_lists.n)+x_lists.m(x_lists.n));
        
        J = G1 + G2;
        J = bsxfun(@times, J(x_lists.n,:), sign_factor) + G3(x_lists.p, :);
        J = [J, 2*lambda*G_m2(x_lists.p,:)];
        J(x_lists.m(x_lists.p) > 0,:) = 2*J(x_lists.m(x_lists.p) > 0, :);
        
        J = [real(J); imag(J)].';
        g = 2*J.'*F;
        
        H_plus = cell(size(H1)); H_minus = cell(size(H1));
        for ii = 1:length(H_plus)
            tmp = (sign_factor*sign_factor.').*(H3{ii}(x_lists.n,x_lists.n) + H3{ii}(x_lists.n,x_lists.n).') ...
                + bsxfun(@times, H1{ii}(x_lists.p, x_lists.n), sign_factor.') + bsxfun(@times, H1{ii}(x_lists.p, x_lists.n).', sign_factor)...
                + bsxfun(@times, H2{ii}(x_lists.p, x_lists.n), sign_factor.') + bsxfun(@times, H2{ii}(x_lists.p, x_lists.n).', sign_factor);
            
            tmp2 = bsxfun(@times, H3{ii}(x_lists.p, x_lists.n), sign_factor.') + bsxfun(@times, H3{ii}(x_lists.n, x_lists.p), sign_factor).' ...
                + (sign_factor*sign_factor.').*H1{ii}(x_lists.n,x_lists.n) + H1{ii}(x_lists.p, x_lists.p).'...
                + (sign_factor*sign_factor.').*H2{ii}(x_lists.n,x_lists.n) + H2{ii}(x_lists.p, x_lists.p).';
            
            H_plus{ii} = tmp + tmp2;
            H_minus{ii} = tmp - tmp2;
            
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
        m_zero_inds = x_lists.m(x_lists.p)==0;
        H_m2_real = H_m2; H_m2_imag = H_m2;
        for ii = 1:length(H_m2)
            H_m2_real{ii}(m_zero_inds, m_zero_inds) = (H_m2_real{ii}(m_zero_inds, m_zero_inds) + ...
                bsxfun(@times, H_m2_real{ii}(m_zero_inds, m_zero_inds), sign_factor(m_zero_inds)))/2;
            H_m2_imag{ii}(m_zero_inds, m_zero_inds) = (H_m2_imag{ii}(m_zero_inds, m_zero_inds) - ...
                bsxfun(@times, H_m2_imag{ii}(m_zero_inds, m_zero_inds), sign_factor(m_zero_inds)))/2;
        end
        H_m2 = cellfun(@(x,y) 2*lambda*[real(x), zeros(size(x)); zeros(size(x)), real(y)], H_m2_real, H_m2_imag, 'UniformOutput', 0);
        H = [H(:);H_m2(:)];
        
        h = 0;
        for ii = 1:length(H)
            h = h + F(ii)*H{ii};
        end
        h = h + J.'*J;
        h = 2*h;
end