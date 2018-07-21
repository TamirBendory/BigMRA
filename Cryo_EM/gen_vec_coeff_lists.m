function [L_list, m_list, s_list, p_inds, n_inds] = gen_vec_coeff_lists(maxL, s_lens)

L_list = []; m_list = []; s_list = [];
for l = 0:maxL
    L_list = [L_list; l*ones(s_lens(l+1)*(2*l+1), 1)];
    m_mat = repmat(-l:l, s_lens(l+1), 1);
    m_list = [m_list; m_mat(:)];
    s_mat = repmat((1:s_lens(l+1)).', 1, 2*l+1);
    s_list = [s_list; s_mat(:)];
end

p_inds = find(m_list >= 0);
n_inds = zeros(size(p_inds));
for ii = 1:length(p_inds)
    n_inds(ii) = find(L_list == L_list(p_inds(ii)) & m_list == -m_list(p_inds(ii)) & s_list == s_list(p_inds(ii)));
end