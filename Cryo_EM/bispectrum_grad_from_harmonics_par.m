function [G1, G2, G3] = bispectrum_grad_from_harmonics_par(a_lms, L, B, W, L1_list, L2_list, L3_list, L_list, m_list, s_lens)


if ~exist('W', 'var') || isempty(W)
    W = precomp_wigner_weights(length(a_lms)-1);
end

if ~exist('B','var') || isempty(B)
    B = precomp_bispectrum_coeffs(length(a_lms)-1, size(a_lms{1},1), L, W);
    [B, L1_list, L2_list, L3_list] = form_B_cell(B);
end

maxK = length(B)-1;
G1 = cellfun(@(y) zeros(length(L_list), size(y{1},1), size(y{1},1)), B, 'UniformOutput', 0);
G2 = G1; G3 = G1;
for k = 0:maxK
    B_k = B{k+1};
    q_k = size(B_k{1},1);
    G1_acc = 0; G2_acc = 0; G3_acc = 0;
    parfor ii = 1:length(B_k)
        L1 = L1_list(ii); s1_curr = s_lens(L1+1); 
        L2 = L2_list(ii); s2_curr = s_lens(L2+1);
        jj_3 = L3_list(ii); L3 = abs(L1-L2) + jj_3 - 1; s3_curr = s_lens(L3+1);
        
        B_k_curr = B_k{ii}(:,:, 1:s1_curr, 1:s2_curr, 1:s3_curr);
        % Perform summation over s1, s2, s3:
        tmp_1 = conj(reshape(B_k_curr, q_k^2*s1_curr*s2_curr, s3_curr))*a_lms{L3+1}; % (q1, q2, s1, s2) x m3
        tmp_1 = reshape(tmp_1.', [2*L3+1, q_k, q_k, s1_curr, s2_curr]); % m3 x q1 x q2 x s1 x s2
        tmp_2 = permute(tmp_1, [1,5,2,3,4]); % m3 x s2 x q1 x q2 x s1
        
        tmp_3 = permute(B_k_curr, [5,1,2,3,4]); % s3 x q1 x q2 x s1 x s2
        
        tmp_1 = reshape(tmp_1, (2*L3+1)*q_k^2*s1_curr, s2_curr)*conj(a_lms{L2+1}); % (m3 x q1 x q2 x s1) x m2
        tmp_1 = reshape(tmp_1.', [2*L2+1, 2*L3+1, q_k, q_k, s1_curr]); % m2 x m3 x q1 x q2 x s1
        tmp_1 = permute(tmp_1, [1,2,5,3,4]); % m2 x m3 x s1 x q1 x q2
        
        tmp_3 = reshape(tmp_3, s3_curr*q_k^2*s1_curr, s2_curr)*a_lms{L2+1}; % s3 x q1 x q2 x s1 x m2
        tmp_3 = reshape(tmp_3.', [2*L2+1, s3_curr, q_k, q_k, s1_curr]); % m2 x s3 x q1 x q2 x s1
        
        tmp_2 = reshape(tmp_2, (2*L3+1)*s2_curr*q_k^2, s1_curr)*conj(a_lms{L1+1}); % (m3 x s2 x q1 x q2) x m1
        tmp_2 = reshape(tmp_2.', [2*L1+1, 2*L3+1, s2_curr, q_k, q_k]); % m1 x m3 x s2 x q1 x q2
        
        tmp_3 = reshape(tmp_3, (2*L2+1)*s3_curr*q_k^2, s1_curr)*a_lms{L1+1}; % (m2, s3, q1, q2) x m1
        tmp_3 = reshape(tmp_3.', [2*L1+1, 2*L2+1, s3_curr, q_k, q_k]); % m1 x m2 x s3 x q1 x q2
        
        % Perform summation over m2, m3:
        g1 = zeros(length(L_list), q_k, q_k);
        g2 = zeros(length(L_list), q_k, q_k);
        g3 = zeros(length(L_list), q_k, q_k);
        for m2 = -L2:L2
            for m1 = max(-L1, -L3-m2):min(L1, L3-m2)
                m3 = m1+m2;
                g1(L_list == L1 & m_list == m1, :, :) = g1(L_list == L1 & m_list == m1, :, :)...
                    + conj(W{L1+1,L2+1}(jj_3, m1+L1+1, m2+L2+1))*reshape(tmp_1(m2+L2+1, m3+L3+1, :,:,:), s1_curr, q_k, q_k); % s1 x q1 x q2
                
                g2(L_list == L2 & m_list == m2, :, :) = g2(L_list == L2 & m_list == m2, :, :)...
                    + conj(W{L1+1,L2+1}(jj_3, m1+L1+1, m2+L2+1))*reshape(tmp_2(m1+L1+1, m3+L3+1, :,:,:), s2_curr, q_k, q_k); % s2 x q1 x q2
                
                g3(L_list == L3 & m_list == m3, :, :) = g3(L_list == L3 & m_list == m3, :, :)...
                    + W{L1+1,L2+1}(jj_3, m1+L1+1, m2+L2+1)*reshape(tmp_3(m1+L1+1, m2+L2+1, :,:,:), s3_curr, q_k, q_k); % s3 x q1 x q2
            end
        end
        G1_acc = G1_acc + g1;
        G2_acc = G2_acc + g2;
        G3_acc = G3_acc + g3;
    end
    G1{k+1} = G1_acc;
    G2{k+1} = G2_acc;
    G3{k+1} = G3_acc;
end
end
