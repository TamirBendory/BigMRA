function [x_aligned,max_corr, ind_out] = align_to_reference(x, xref)

    assert(all(size(x) == size(xref)), 'x and xref must have identical size');
    x = x(:);
    xref = xref(:);

    x_fft = fft(x);
    xref_fft = fft(xref);
    
    correlation_x_xref = real(ifft(conj(x_fft) .* xref_fft));
 
    [max_corr, ind] = max(correlation_x_xref);
    
    ind_out = ind - 1;
    
    max_corr = max_corr/norm(xref_fft);
 
    x_aligned = circshift(x, ind-1);

end
