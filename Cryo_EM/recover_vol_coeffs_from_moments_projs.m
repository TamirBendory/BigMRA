function a_lms = recover_vol_coeffs_from_moments_projs(m3_micro, maxL, L, resol, r_cut, gamma)

s_lens = gen_s_list(maxL, r_cut, 1, floor(L/2));
L_list = gen_vec_coeff_lists(maxL, s_lens);
Rots = genRotationsGrid(resol, 0, 0);

problem.M = euclideanfactory(length(L_list), 1);
problem.cost = @(a_lms, store) costgrad_micro_moments_projs(a_lms, L_list, L, Rots, m3_micro, gamma, store);

% checkgradient(problem)

a_lms = conjugategradient(problem, [], []);
a_lms = accumarray(L_list, a_lms, [max(L_list)+1, 1], @(x) {x});
for l = 0:length(a_lms)-1
    a_lms{l+1} = reshape(a_lms{l+1}, [], 2*l+1);
end