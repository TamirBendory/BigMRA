function [m1, m2, m3] = moments_from_micrograph_steerable_windows_GPU(I, L)
% L is size of single projection
m1 = mean(I(:))/L;

beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
T = 1e-1;       % Truncation parameter

[Wt, ang_freq] = precomp_pswf_t_windows(L-1, beta, T);
Wt = gpuArray(Wt);

maxN = max(ang_freq);
q_list = zeros(maxN+1, 1);
for ii = 0:maxN
    q_list(ii+1) = sum(ang_freq == ii); 
end
q_cumsum = cumsum([0; q_list(:)]);
q_sq_cumsum = cumsum([0; q_list(:).^2]);

m2 = zeros(q_list(1), 1, 'double');

blk_id = [];
for ii = 0:maxN
    blk_id(end+1:end+q_list(ii+1)^2, 1) = ii+1;
end
m3 = zeros(size(blk_id), 'double');

sz_img = size(I);
pixels_per_img = size(I,1)*size(I,2); % pixels per micrograph
num_micros = size(I,3);

for t = 1:num_micros
    I_curr = gpuArray(I(:,:,t));
    I_pad = fft2(I_curr, sz_img(1)+2*L-2, sz_img(2)+2*L-2);
    
    m3_add = zeros(length(blk_id), 1, 'double');
    for N = 0:maxN
        coeff_N = zeros(sz_img(1), sz_img(2), q_list(N+1), 'double');
        
        for ii = q_cumsum(N+1)+1:q_cumsum(N+2)
            coeff = ifft2(I_pad.*fft2(Wt(:,:,ii), sz_img(1)+2*L-2, sz_img(2)+2*L-2));
            coeff = gather(coeff);
%             coeff = xcorr2_fft(I_curr, Wt(:,:, ii)); % expansion in PSWFs
            coeff_N(:,:,ii-q_cumsum(N+1)) = coeff(L:end-L+1, L:end-L+1);
        end
        
        coeff_N = reshape(coeff_N, pixels_per_img, q_list(N+1)).';
        if N == 0
            m2 = m2 + real(coeff_N*I_curr(:));
        end
            
        tmp = bsxfun(@times, coeff_N, I_curr(:).');
        tmp = real(tmp*coeff_N');
        m3_add(q_sq_cumsum(N+1)+1: q_sq_cumsum(N+2)) = tmp(:);
    end
    m3 = m3 + m3_add;
end
    

m2 = m2./(pixels_per_img*num_micros);
m3 = m3./(pixels_per_img*num_micros);

m3 = accumarray(blk_id, m3, [max(ang_freq)+1, 1], @(x) {reshape(x, sqrt(length(x)), [])});
