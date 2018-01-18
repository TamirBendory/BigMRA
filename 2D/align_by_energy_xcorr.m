function x_out = align_by_energy_xcorr(x,L)
% input: 
% x  -  the zero-padded input signal
% L  -  size of output image
% output: 
% x_out - L X L subimage of x that contains most of the energy

I = xcorr2(x.^2,ones(L));
[~,ind] = max(I(:));
[ind1, ind2] = ind2sub(size(I),ind);
x_out = x(ind1-L+1:ind1,ind2-L+1:ind2);