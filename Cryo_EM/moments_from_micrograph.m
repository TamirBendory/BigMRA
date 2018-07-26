function [m1, m2, m3] = moments_from_micrograph(I, L)

m1 = mean(I(:))/L;

m2 = zeros(2*L-1, 2*L-1);
parfor ii = 1:(2*L-1)^2
    [di, dj] = ind2sub([2*L-1, 2*L-1], ii);
    di = di - L; dj = dj - L;
    I_prod = I .* aperiodic_shift(I, [di, dj]);
    m2(ii) = sum(I_prod(:))./numel(I);
end

m3 = zeros(2*L-1, 2*L-1, 2*L-1, 2*L-1);
parfor ii = 1:(2*L-1)^4
    [di1, dj1, di2, dj2] = ind2sub([2*L-1, 2*L-1, 2*L-1, 2*L-1], ii);
    di1 = di1 - L; dj1 = dj1 - L; di2 = di2 - L; dj2 = dj2 - L;
    I_prod = I .* aperiodic_shift(I, [di1, dj1]) .* aperiodic_shift(I, [di2, dj2]);
    m3(ii) = sum(I_prod(:))./numel(I);
end
    
% m2 = zeros(2*L-1, 2*L-1); m3 = m2;
% I_spec = aperiodic_shift(I, [1, 1]);
% parfor ii = 1:(2*L-1)^2
%     [di, dj] = ind2sub([2*L-1, 2*L-1], ii);
%     di = di - L; dj = dj - L;
%     I_prod = I .* aperiodic_shift(I, [di, dj]);
%     m2(ii) = sum(I_prod(:))./numel(I);
%     m3(ii) = sum(sum(sum(I_prod.*I_spec)))./numel(I);
% end
