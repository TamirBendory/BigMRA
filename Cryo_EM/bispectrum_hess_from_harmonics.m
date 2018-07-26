function [H1, H2, H3] = bispectrum_hess_from_harmonics(a_lms, L, B, W, L1_list, L2_list, L3_list, L_list, m_list, s_lens, blk_id)

if ~exist('W', 'var') || isempty(W)
    W = precomp_wigner_weights(length(a_lms)-1);
end

if ~exist('B','var') || isempty(B)
    [B, L1_list, L2_list, L3_list, blk_id] = precomp_bispectrum_coeffs(length(a_lms)-1, s_lens, L, W);
end

kq_len = length(blk_id);
H1(length(L_list), length(L_list), kq_len) = 0; 
H2(length(L_list), length(L_list), kq_len) = 0; 
H3(length(L_list), length(L_list), kq_len) = 0;
parfor ii = 1:length(B)
    L1 = L1_list(ii); s1_curr = s_lens(L1+1);
    L2 = L2_list(ii); s2_curr = s_lens(L2+1);
    jj_3 = L3_list(ii); L3 = abs(L1-L2) + jj_3 - 1; s3_curr = s_lens(L3+1);
    
    % Perform summation over s1, s2, s3 (re-using tmp_1 to decrease memory and allocation costs):
    tmp_1 = reshape(B{ii}(:, 1:s3_curr, 1:s2_curr, 1:s1_curr), kq_len*s3_curr*s2_curr, s1_curr)*a_lms{L1+1}; % k x s3 x s2 x m1
    tmp_2 = reshape(permute(B{ii}(:, 1:s3_curr, 1:s2_curr, 1:s1_curr), [1,2,4,3]), kq_len*s3_curr*s1_curr, s2_curr)*a_lms{L2+1}; % k x s3 x s1 x m2
    tmp_3 = reshape(permute(B{ii}(:, 1:s3_curr, 1:s2_curr, 1:s1_curr), [1,3,4,2]), kq_len*s2_curr*s1_curr, s3_curr)*conj(a_lms{L3+1}); % k x s2 x s1 x m3
    
    h1 = zeros(length(L_list), length(L_list), kq_len, 'double'); 
    h2 = zeros(length(L_list), length(L_list), kq_len, 'double');  
    h3 = zeros(length(L_list), length(L_list), kq_len, 'double'); 
    for m1 = -L1:L1
        for m2 = max(-L2, -L3-m1):min(L2, L3-m1)
            m3 = m1+m2;
            norm_fact = W{L1+1,L2+1}(m2+L2+1, m1+L1+1, jj_3);
            
            h1(L_list == L3 & m_list == m3, L_list == L2 & m_list == m2, :) = h1(L_list == L3 & m_list == m3, L_list == L2 & m_list == m2, :) + ...
                norm_fact*permute(reshape(tmp_1(:, m1+L1+1), kq_len, s3_curr, s2_curr), [2,3,1]); % s3 x s2 x k
            
            h2(L_list == L3 & m_list == m3, L_list == L1 & m_list == m1, :) = h2(L_list == L3 & m_list == m3, L_list == L1 & m_list == m1, :) + ...
                norm_fact*permute(reshape(tmp_2(:, m2+L2+1), kq_len, s3_curr, s1_curr), [2,3,1]); % s3 x s1 x k
            
            h3(L_list == L2 & m_list == m2, L_list == L1 & m_list == m1, :) = h3(L_list == L2 & m_list == m2, L_list == L1 & m_list == m1, :) + ...
                norm_fact*permute(reshape(tmp_3(:, m3+L3+1), kq_len, s2_curr, s1_curr), [2,3,1]); % s2 x s1 x k
        end
    end
    H1 = H1 + h1;
    H2 = H2 + h2;
    H3 = H3 + h3;
end
H1 = mat2cell(H1, length(L_list), length(L_list), ones(kq_len, 1));
H2 = mat2cell(H2, length(L_list), length(L_list), ones(kq_len, 1));
H3 = mat2cell(H3, length(L_list), length(L_list), ones(kq_len, 1));
