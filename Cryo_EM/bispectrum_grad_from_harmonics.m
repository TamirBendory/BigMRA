function [G1, G2, G3] = bispectrum_grad_from_harmonics(a_lms, L, B, W)

maxL = length(a_lms)-1;
maxS = size(a_lms{1},1);

if ~exist('W', 'var') || isempty(W)
    W = precomp_wigner_weights(maxL);
end

if ~exist('B','var') || isempty(B)
    B = precomp_bispectrum_coeffs(maxL, maxS, L, W);
end

s_lens = cellfun('size', a_lms, 1);
[L_list, m_list, ~] = gen_vec_coeff_lists(maxL, s_lens);

q_list = cellfun(@(x) size(x,1), B{1}{1}, 'UniformOutput', 0);
G1 = cellfun(@(y) zeros(length(L_list), y, y), q_list, 'UniformOutput', 0);
G2 = G1; G3 = G1;
for ii = 1:(maxL+1)^2
    [L1, L2] = ind2sub([maxL+1, maxL+1], ii);
    L1 = L1-1; L2 = L2-1;
    
    L3_vals = abs(L1-L2):min(L1+L2, maxL);
    for jj = 1:length(L3_vals)
        L3 = L3_vals(jj);
        
        for k = 0:length(B{ii}{jj})-1
            sz_B = size(B{ii}{jj}{k+1}); % q1 x q2 x s1 x s2 x s3
            % Perform summation over s1, s2, s3:
            sz = size(a_lms{L3+1});
            tmp = conj(reshape(B{ii}{jj}{k+1}(:,:,:,:,1:sz(1)), [], sz(1)))*a_lms{L3+1}; % (q1, q2, s1, s2) x m3
            tmp = reshape(tmp.', [sz(2), sz_B(1:4)]); % m3 x q1 x q2 x s1 x s2
            tmp_3 = permute(B{ii}{jj}{k+1}(:,:,:,:,1:s_lens(L3+1)), [5,1,2,3,4]); % s3 x q1 x q2 x s1 x s2
            
            tmp_1 = tmp;
            tmp_2 = permute(tmp(:,:,:,:,1:s_lens(L2+1)), [1,5,2,3,4]); % m3 x s2 x q1 x q2 x s1
            
            
            sz_B = size(tmp_1);
            sz = size(a_lms{L2+1});
            tmp_1 = reshape(tmp_1(:,:,:,:,1:sz(1)), [], sz(1))*conj(a_lms{L2+1}); % (m3 x q1 x q2 x s1) x m2
            tmp_1 = reshape(tmp_1.', [sz(2), sz_B(1:4)]); % m2 x m3 x q1 x q2 x s1
            tmp_1 = permute(tmp_1(:,:,:,:,1:s_lens(L1+1)), [1,2,5,3,4]); % m2 x m3 x s1 x q1 x q2
            
            sz_G = size(tmp_3);
            tmp_3 = reshape(tmp_3(:,:,:,:,1:sz(1)), [], sz(1))*a_lms{L2+1}; % s3 x q1 x q2 x s1 x m2
            tmp_3 = reshape(tmp_3.', [sz(2), sz_G(1:4)]); % m2 x s3 x q1 x q2 x s1
            
            sz_B = size(tmp_2);
            sz = size(a_lms{L1+1});
            tmp_2 = reshape(tmp_2(:,:,:,:,1:sz(1)), [], sz(1))*conj(a_lms{L1+1}); % (m3 x s2 x q1 x q2) x m1
            tmp_2 = reshape(tmp_2.', [sz(2), sz_B(1:4)]); % m1 x m3 x s2 x q1 x q2
            
            sz_G = size(tmp_3);
            sz = size(a_lms{L1+1});
            tmp_3 = reshape(tmp_3(:,:,:,:,1:sz(1)), [], sz(1))*a_lms{L1+1}; % (m2, m3, q1, q2) x m1
            tmp_3 = reshape(tmp_3.', [sz(2), sz_G(1:4)]); % m1 x m2 x s3 x q1 x q2
            
            % Perform summation over m2, m3:
            for m1 = -L1:L1
                for m2 = -L2:L2
                    m3 = m1+m2;
                    if abs(m3) > L3, continue; end
                    G1{k+1}(L_list == L1 & m_list == m1, :, :) = G1{k+1}(L_list == L1 & m_list == m1, :, :)...
                        + conj(W{ii}(jj, m1+L1+1, m2+L2+1))*squeeze(tmp_1(m2+L2+1, m3+L3+1, :,:,:)); % s1 x q1 x q2
                    G2{k+1}(L_list == L2 & m_list == m2, :, :) = G2{k+1}(L_list == L2 & m_list == m2, :, :)...
                        + conj(W{ii}(jj, m1+L1+1, m2+L2+1))*squeeze(tmp_2(m1+L1+1, m3+L3+1, :,:,:)); % s2 x q1 x q2
                    G3{k+1}(L_list == L3 & m_list == m3, :, :) = G3{k+1}(L_list == L3 & m_list == m3, :, :)...
                        + W{ii}(jj, m1+L1+1, m2+L2+1)*squeeze(tmp_3(m1+L1+1, m2+L2+1, :,:,:)); % s3 x q1 x q2
                end
            end
        end
    end
end
