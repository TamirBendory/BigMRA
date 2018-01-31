function X_aligned = align_to_reference_2D(X, X_ref)

    C = ifft2(conj(fft2(X_ref)) .* fft2(X));
    [~, ind] = max(C(:));
    [k1, k2] = ind2sub(size(C), ind);
    X_aligned = circshift(X, [1-k1, 1-k2]);

end
