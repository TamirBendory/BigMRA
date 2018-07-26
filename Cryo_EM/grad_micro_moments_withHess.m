function [g, store] = grad_micro_moments_withHess(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda, store)

if ~isfield(store, 'F')
    [~, store] = cost_micro_moments_withHess(a_lms, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda, store);
else
    G1 = store.G1;
    G2 = store.G2;
    G3 = store.G3;
    G_m2 = store.G_m2;
    F = store.F;
end

if ~isfield(store, 'J')
    sign_factor = (-1).^(x_lists.L(x_lists.n)+x_lists.m(x_lists.n));
    
    J = G1 + G2; 
    J = bsxfun(@times, J(x_lists.n,:), sign_factor) + G3(x_lists.p, :);
    J = [J, 2*lambda*G_m2(x_lists.p,:)];
    J(x_lists.m(x_lists.p) > 0,:) = 2*J(x_lists.m(x_lists.p) > 0, :);
    
    J = [real(J); imag(J)].';
    
    store.J = J;
    store.signs = sign_factor;
else
    J = store.J;
end

g = 2*J.'*F;