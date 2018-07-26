function psi_Nn = PSWF_2D_full_cart_OLD(maxN, n_list, L)

[x,y] = meshgrid(-L:L-1, -L:L-1);
[phi,r] = cart2pol(x,y);
r = r./L;

psi_Nn = cell(2*maxN+1, 1);
for N = -maxN:maxN
    psi_Nn{N+maxN+1} = zeros((2*L)^2, n_list(abs(N)+1));
    R_Nn = PSWF_radial_2D(abs(N), n_list(abs(N)+1)-1, pi*L, r(r<=1)); % generate 2D radial prolates
    psi_Nn{N+maxN+1}(r<=1,:) = bsxfun(@times, R_Nn, exp(1i*N*phi(r<=1))./sqrt(2*pi));
end
    
    