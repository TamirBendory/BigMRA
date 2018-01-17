function x = circshift_ad(x, k)
% Quick implementation of circshift for vectors, because Matlab's circshift
% is not recognized by ADiMat (for automatic differentiation in Manopt.)

    n = size(x, 1);
    k = mod(k, n);
    x = x([end-(k-1):end, 1:(end-k)], :);

end
