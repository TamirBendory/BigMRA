function x_out = align_by_energy(x, L)
% input: 
%    x  -  an image (matrix) of size WxW, with W >= L
%    L  -  desired output size
%
% output: 
%    x_out - LxL subimage of x that contains most of the energy (in Frobenius norm)
    
    I = conv2(x.*conj(x), ones(L))
    [~, ind] = max(I(:));
    [ind1, ind2] = ind2sub(size(I), ind);
    x_out = x(ind1-L+1:ind1, ind2-L+1:ind2);
    
end
