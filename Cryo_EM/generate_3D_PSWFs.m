function [Psi_lms, pts_in_ball, Psi_lms_2D, pts_in_disc_2D] = generate_3D_PSWFs(L, maxL)

R = floor(L/2);
c=pi*R;   % band limit
D=3;      % prolate on 3-D ball
minEigenvalRatio = 10^-16; % the truncation parameter. magnitude of the smallest eigenvalue that we would like to keep
matdim = 1000;                        % size of the matrix representation of the differential operator, see notes. If this is too small, a warning should appear. This will be removed in future versions.    
prolate_crea_options.isfixfirst = 1; % this option should be set to 1 for more accurate eigenvalues.

if mod(L,2) == 0
    [x,y,z] = meshgrid(-R:R-1, -R:R-1, -R:R-1);
else
    [x,y,z] = meshgrid(-R:R, -R:R, -R:R);
end
r = sqrt( x.^2 + y.^2 + z.^2 );
theta = acos( z./r);
phi = acos( x./(r.*sin(theta)) );
j1 = find( y < 0);
phi(j1) = 2*pi - phi(j1) ;
MACHINEEPS = 1e-15;
j = find( theta < MACHINEEPS | theta > pi-MACHINEEPS); % z = 1 or -1, x=y=0
phi(j) = 0;
phi = real( phi); %enforce real numbers
jorigin = find( r < MACHINEEPS ); % be careful at orgin point
theta(jorigin) = 0;
phi( jorigin ) = 0;
r = r/R;

pts_in_ball = find(r <= 1);
pts_in_disc_2D = find(r(z==0) <= 1);

Psi_lms = cell(maxL+1, 1);
Psi_lms_2D = cell(maxL+1, 1);
Y_l = YN2YL(getSH(maxL, [phi(pts_in_ball), theta(pts_in_ball)], 'complex'));
for l = 0:maxL
    prolate_dat = prolate_crea(c,D,l,minEigenvalRatio, matdim, prolate_crea_options);
    Phi_ls = prolate_ev(prolate_dat, 0:prolate_dat.num_prols-1 , r(pts_in_ball));
    Psi_lms{l+1} = bsxfun(@times, Phi_ls, permute(Y_l{l+1},[2,3,1]));
    Psi_lms_2D{l+1} = Psi_lms{l+1}(z(pts_in_ball) == 0, :, 1:2:end); % restrict to xy-plane and take only m = l (mod 2)
end

end